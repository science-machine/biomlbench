#!/usr/bin/env python3
import argparse
import asyncio
import json
import logging
import os
import random
import signal
import sys
import time
from dataclasses import dataclass
from itertools import cycle
from typing import Iterable, List, Optional

from tenacity import (  # type: ignore
    RetryError,
    retry,
    retry_if_exception_type,
    stop_after_attempt,
    wait_random_exponential,
)

# Lazy import anthropic so the script can print a helpful error if missing
try:
    from anthropic import AsyncAnthropic, APIStatusError, RateLimitError
except Exception:  # pragma: no cover - only for helpful error message
    AsyncAnthropic = None  # type: ignore
    APIStatusError = Exception  # type: ignore
    RateLimitError = Exception  # type: ignore

from dotenv import load_dotenv  # type: ignore

load_dotenv()

@dataclass
class Args:
    model: str
    max_output_tokens: int
    concurrency: int
    min_sleep: float
    max_sleep: float
    places_file: Optional[str]
    use_thinking: bool
    thinking_budget_tokens: int
    temperature: float
    system_prompt: str
    max_requests_per_worker: int
    log_json: bool


DEFAULT_PLACES: List[str] = [
    "Kyoto, Japan",
    "Reykjavik, Iceland",
    "Cusco, Peru",
    "Marrakesh, Morocco",
    "Queenstown, New Zealand",
    "Barcelona, Spain",
    "Banff, Canada",
    "Santorini, Greece",
    "Seoul, South Korea",
    "Cape Town, South Africa",
    "Hoi An, Vietnam",
    "Amalfi Coast, Italy",
    "Istanbul, Türkiye",
    "Petra, Jordan",
    "Bali, Indonesia",
    "Lisbon, Portugal",
    "Tulum, Mexico",
    "Prague, Czech Republic",
    "Edinburgh, Scotland",
    "Patagonia, Argentina",
]


SYSTEM_PROMPT_DEFAULT = (
    "You are an expert travel planner. Produce practical, detailed, and creative itineraries. "
    "Focus on logistics, budgets, transport, and local insights."
)


def build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Spend Anthropic credits via varied travel-planning prompts.")
    p.add_argument("--model", default="claude-sonnet-4-20250514", help="Anthropic model name")
    p.add_argument("--max-output-tokens", type=int, default=4000, help="Max output tokens per call")
    p.add_argument("--concurrency", type=int, default=5, help="Number of concurrent workers")
    p.add_argument("--min-sleep", type=float, default=0.5, help="Min seconds to sleep between calls per worker")
    p.add_argument("--max-sleep", type=float, default=2.0, help="Max seconds to sleep between calls per worker")
    p.add_argument("--places-file", type=str, default=None, help="Optional path to a file with one place per line")
    p.add_argument("--use-thinking", action="store_true", help="Enable thinking mode if the model supports it")
    p.add_argument("--thinking-budget-tokens", type=int, default=2048, help="Budget for thinking tokens, if enabled")
    p.add_argument("--temperature", type=float, default=0.7, help="Sampling temperature")
    p.add_argument(
        "--system-prompt",
        type=str,
        default=SYSTEM_PROMPT_DEFAULT,
        help="System prompt to use for all requests",
    )
    p.add_argument(
        "--max-requests-per-worker",
        type=int,
        default=0,
        help="Optional cap per worker; 0 means run forever",
    )
    p.add_argument("--log-json", action="store_true", help="Emit logs as JSON for easier ingestion")
    return p


def setup_logging(log_json: bool) -> None:
    level = logging.INFO
    handler = logging.StreamHandler(sys.stdout)
    if log_json:
        formatter = logging.Formatter("%(message)s")
    else:
        formatter = logging.Formatter("[%(asctime)s] %(levelname)s %(message)s")
    handler.setFormatter(formatter)
    root = logging.getLogger()
    root.setLevel(level)
    root.handlers.clear()
    root.addHandler(handler)


def load_places(places_file: Optional[str]) -> List[str]:
    if places_file is None:
        return DEFAULT_PLACES
    try:
        with open(places_file, "r", encoding="utf-8") as f:
            places = [line.strip() for line in f if line.strip()]
        return places or DEFAULT_PLACES
    except FileNotFoundError:
        return DEFAULT_PLACES


def random_travel_context() -> dict:
    durations = ["3 days", "5 days", "1 week", "10 days", "2 weeks"]
    budgets = [
        "$800", "$1,500", "$2,500", "$4,000", "$6,000", "$10,000",
        "mid-range", "luxury", "backpacker",
    ]
    styles = [
        "food-focused",
        "adventure",
        "slow travel",
        "family-friendly",
        "nightlife",
        "art and culture",
        "nature and hiking",
        "photography",
        "wellness and spa",
    ]
    seasons = [
        "early spring",
        "late spring",
        "summer",
        "early autumn",
        "late autumn",
        "winter",
    ]
    transport = [
        "public transit",
        "rental car",
        "walkable core",
        "mix of transit and ride-shares",
    ]
    return {
        "duration": random.choice(durations),
        "budget": random.choice(budgets),
        "style": random.choice(styles),
        "season": random.choice(seasons),
        "transport": random.choice(transport),
    }


def build_user_prompt(place: str) -> str:
    ctx = random_travel_context()
    variants = [
        (
            "Plan an in-depth itinerary for 2 people to {place} over {duration} in {season}. "
            "Prioritize {style} experiences with a {budget} budget. Assume {transport}. "
            "Include daily schedules, approximate costs, intra-city transit, and 2-3 backup options."
        ),
        (
            "Create a logistics-first travel plan for a couple visiting {place} for {duration} in {season}. "
            "Focus on routes, neighborhoods to stay, realistic timing, and {style} highlights within a {budget} budget. "
            "Add booking tips, safety notes, and low-cost alternates."
        ),
        (
            "Draft a detailed travel brief for two travelers going to {place}. Trip length: {duration}; season: {season}; style: {style}; budget: {budget}. "
            "Provide a structured plan with morning/afternoon/evening blocks, travel time estimates, and dining suggestions."
        ),
        (
            "Design a day-by-day plan for {duration} in {place} during {season} for 2 people with a {style} focus and {budget} budget. "
            "Add pitfalls to avoid, must-try local experiences, and a packing checklist."
        ),
    ]
    template = random.choice(variants)
    return template.format(place=place, **ctx)


def json_log(enabled: bool, **fields: object) -> None:
    if not enabled:
        logging.info(
            " ".join(f"{k}={v}" for k, v in fields.items())
        )
        return
    try:
        logging.info(json.dumps(fields, ensure_ascii=False))
    except Exception:
        logging.info(str(fields))


def _get_client(args: Args) -> AsyncAnthropic:
    if AsyncAnthropic is None:
        print(
            "Missing dependency 'anthropic'. Install it with: pip install anthropic",
            file=sys.stderr,
        )
        sys.exit(2)
    api_key = os.getenv("ANTHROPIC_API_KEY")
    if not api_key:
        print(
            "Environment variable ANTHROPIC_API_KEY is not set. Put it in your .env or environment.",
            file=sys.stderr,
        )
        sys.exit(2)
    default_headers = {"anthropic-beta": "thinking-v1"} if args.use_thinking else None
    return AsyncAnthropic(api_key=api_key, timeout=60, default_headers=default_headers)


class AnthropicSender:
    def __init__(self, client: AsyncAnthropic, args: Args) -> None:
        self._client = client
        self._args = args

    @retry(
        reraise=True,
        stop=stop_after_attempt(5),
        wait=wait_random_exponential(multiplier=1, max=30),
        retry=retry_if_exception_type((RateLimitError, APIStatusError, asyncio.TimeoutError)),
    )
    async def send_once(self, place: str) -> dict:
        user_prompt = build_user_prompt(place)
        request_kwargs = {
            "model": self._args.model,
            "max_tokens": self._args.max_output_tokens,
            "temperature": self._args.temperature,
            "system": self._args.system_prompt,
            "messages": [
                {"role": "user", "content": user_prompt},
            ],
        }
        # Optional thinking support (best-effort)
        if self._args.use_thinking:
            request_kwargs["extra_body"] = {
                "thinking": {
                    "type": "enabled",
                    "budget_tokens": self._args.thinking_budget_tokens,
                }
            }
        try:
            resp = await self._client.messages.create(**request_kwargs)  # type: ignore[arg-type]
            # Extract token usage if available
            usage = getattr(resp, "usage", None)
            input_tokens = getattr(usage, "input_tokens", None) if usage else None
            output_tokens = getattr(usage, "output_tokens", None) if usage else None
            return {
                "place": place,
                "input_tokens": input_tokens,
                "output_tokens": output_tokens,
                "id": getattr(resp, "id", None),
            }
        except Exception:
            # Re-raise to trigger tenacity retry where applicable
            raise


async def worker(worker_id: int, sender: AnthropicSender, places_iter: Iterable[str], args: Args, stop_event: asyncio.Event) -> None:
    rng = random.Random(time.time() + worker_id)
    count = 0
    for place in places_iter:
        if stop_event.is_set():
            return
        t0 = time.time()
        try:
            result = await sender.send_once(place)
            dt = time.time() - t0
            json_log(
                args.log_json,
                event="request_ok",
                worker=worker_id,
                place=place,
                ms=round(dt * 1000),
                input_tokens=result.get("input_tokens"),
                output_tokens=result.get("output_tokens"),
                id=result.get("id"),
            )
        except RetryError as e:
            # Retries exhausted
            dt = time.time() - t0
            last_exc = e.last_attempt.exception() if e.last_attempt else None
            json_log(
                args.log_json,
                event="request_failed_retried_out",
                worker=worker_id,
                place=place,
                ms=round(dt * 1000),
                error=str(last_exc) if last_exc else "unknown",
            )
        except Exception as e:
            dt = time.time() - t0
            json_log(
                args.log_json,
                event="request_failed",
                worker=worker_id,
                place=place,
                ms=round(dt * 1000),
                error=str(e),
            )
        # Per-worker pacing with jitter
        sleep_s = rng.uniform(args.min_sleep, args.max_sleep)
        await asyncio.sleep(sleep_s)
        count += 1
        if args.max_requests_per_worker > 0 and count >= args.max_requests_per_worker:
            return


async def run(args: Args) -> None:
    setup_logging(args.log_json)

    # Prepare places iterator (cycled to be infinite)
    places = load_places(args.places_file)
    places_iter = cycle(places)

    client = _get_client(args)
    sender = AnthropicSender(client, args)

    stop_event = asyncio.Event()

    def _handle_sigterm() -> None:
        logging.warning("Received termination signal; stopping workers...")
        stop_event.set()

    loop = asyncio.get_running_loop()
    for sig in (signal.SIGINT, signal.SIGTERM):
        try:
            loop.add_signal_handler(sig, _handle_sigterm)
        except NotImplementedError:
            # Windows
            signal.signal(sig, lambda *_: _handle_sigterm())

    tasks = [
        asyncio.create_task(worker(i, sender, places_iter, args, stop_event))
        for i in range(args.concurrency)
    ]

    json_log(
        args.log_json,
        event="runner_started",
        model=args.model,
        concurrency=args.concurrency,
        min_sleep=args.min_sleep,
        max_sleep=args.max_sleep,
        max_output_tokens=args.max_output_tokens,
        use_thinking=args.use_thinking,
        thinking_budget_tokens=args.thinking_budget_tokens,
    )

    await asyncio.gather(*tasks, return_exceptions=True)


def main() -> None:
    parser = build_arg_parser()
    ns = parser.parse_args()
    args = Args(
        model=ns.model,
        max_output_tokens=ns.max_output_tokens,
        concurrency=ns.concurrency,
        min_sleep=ns.min_sleep,
        max_sleep=ns.max_sleep,
        places_file=ns.places_file,
        use_thinking=ns.use_thinking,
        thinking_budget_tokens=ns.thinking_budget_tokens,
        temperature=ns.temperature,
        system_prompt=ns.system_prompt,
        max_requests_per_worker=ns.max_requests_per_worker,
        log_json=ns.log_json,
    )
    if args.min_sleep < 0 or args.max_sleep < 0 or args.max_sleep < args.min_sleep:
        print("Invalid sleep window; ensure 0 <= min <= max", file=sys.stderr)
        sys.exit(2)
    try:
        asyncio.run(run(args))
    except KeyboardInterrupt:
        pass


if __name__ == "__main__":
    main()
