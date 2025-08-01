import os
import re
from pathlib import Path
from typing import Optional


def get_env_var(value: str) -> Optional[str]:
    """Returns the name of the environment variable in the format `${secrets.<name>}`."""

    if not isinstance(value, str):
        return None

    env_var_pattern = r"\$\{\{\s*secrets\.(\w+)\s*\}\}"
    match = re.match(env_var_pattern, value)

    if not match:
        return None

    return match.group(1)


def is_env_var(value: str) -> bool:
    """Checks if the value is an environment variable."""

    return get_env_var(value) is not None


def load_env_file_if_exists() -> bool:
    """
    Loads .env file from the project root if it exists.
    Returns True if file was loaded, False otherwise.
    """
    try:
        from dotenv import load_dotenv
        
        # Find the project root by looking for biomlbench directory
        current_path = Path(__file__).parent
        while current_path != current_path.parent:
            if (current_path / "biomlbench").is_dir():
                env_file = current_path / ".env"
                if env_file.exists():
                    load_dotenv(env_file)
                    return True
                break
            current_path = current_path.parent
        return False
    except ImportError:
        # python-dotenv not available, skip loading
        return False


def parse_env_var_values(dictionary: dict) -> dict:
    """
    Parses any values in the dictionary that match the ${{ secrets.ENV_VAR }} pattern and replaces
    them with the value of the ENV_VAR environment variable.
    
    Provides helpful error messages when environment variables are missing.
    """
    # Try to load .env file if it exists
    load_env_file_if_exists()
    
    for key, value in dictionary.items():
        if not is_env_var(value):
            continue

        env_var = get_env_var(value)
        
        if env_var is None:
            continue  # Skip if not a valid env var pattern

        if os.getenv(env_var) is None:
            # Provide helpful error message with suggestions
            error_msg = f"Environment variable `{env_var}` is not set!"
            
            # Add specific suggestions based on the variable name
            if env_var.endswith("_API_KEY"):
                error_msg += f"\n\nTo fix this:"
                error_msg += f"\n1. Create a .env file in your project root directory"
                error_msg += f"\n2. Add the following line to your .env file:"
                error_msg += f"\n   {env_var}=your-actual-api-key-here"
                error_msg += f"\n3. Or export the variable in your shell:"
                error_msg += f"\n   export {env_var}=your-actual-api-key-here"
            
            raise ValueError(error_msg)

        dictionary[key] = os.getenv(env_var)

    return dictionary
