import logging
import os
from pathlib import Path


ITALIC = "\033[3m"
BOLD = "\033[1m"
RESET = "\033[0m"


def bold(text):
    """Return text in bold."""
    return f"{BOLD}{text}{RESET}"
def italic(text):
    """Return text in italic."""
    return f"{ITALIC}{text}{RESET}"    


class Logger(logging.Logger):
    def __init__(self, name, output_dir=".", output_file_name="log.log"):
        assert os.path.exists(output_dir), f"{output_dir} does not exist, create it first"
        super().__init__(name)
        self.setLevel(logging.INFO)
        # Create handlers
        self.c_handler = logging.StreamHandler()
        self.f_handler = logging.FileHandler(os.path.join(output_dir, output_file_name))
        self.c_handler.setLevel(self.level)
        self.f_handler.setLevel(self.level)
        # Create formatters and add them to handlers
        self.c_format = logging.Formatter(f"[{BOLD}%(levelname)s{RESET}] %(message)s")
        self.f_format = logging.Formatter(f"%(asctime)s {ITALIC}%(filename)s{RESET} [{BOLD}%(levelname)s{RESET}] %(message)s")
        self.c_handler.setFormatter(self.c_format)
        self.f_handler.setFormatter(self.f_format)
        # Add handlers to the logger
        self.addHandler(self.c_handler)
        self.addHandler(self.f_handler)

    def set_debug(self):
        """Set the logger to debug mode."""
        self.setLevel(logging.DEBUG)
        for handler in self.handlers:
            handler.setLevel(logging.DEBUG)
        
    def set_info(self):
        """Set the logger to info mode."""
        self.setLevel(logging.INFO)
        for handler in self.handlers:
            handler.setLevel(logging.INFO)

    def set_warning(self):
        """Set the logger to warning mode."""
        self.setLevel(logging.WARNING)
        for handler in self.handlers:
            handler.setLevel(logging.WARNING)

    def set_error(self):
        """Set the logger to error mode."""
        self.setLevel(logging.ERROR)
        for handler in self.handlers:
            handler.setLevel(logging.ERROR)

    def set_critical(self):
        """Set the logger to critical mode."""
        self.setLevel(logging.CRITICAL)
        for handler in self.handlers:
            handler.setLevel(logging.CRITICAL)

    def change_output_dir_file(self, output_dir, output_file):
        new_handler = logging.FileHandler(os.path.join(output_dir, output_file))
        new_handler.setLevel(self.f_handler.level)
        new_handler.setFormatter(self.f_format)
        self.removeHandler(self.f_handler)
        self.f_handler = new_handler
        self.addHandler(self.f_handler)

    def disable_console(self):
        """Disable the console output."""
        self.removeHandler(self.c_handler)
        self.c_handler.close()
        
    def enable_console(self):
        """Enable the console output."""
        self.addHandler(self.c_handler)
        self.c_handler.setLevel(self.level)
        self.c_handler.setFormatter(self.c_format)

    def disable_file(self):
        """Disable the file output."""
        self.removeHandler(self.f_handler)
        self.f_handler.close()

    def enable_file(self):
        """Enable the file output."""
        self.addHandler(self.f_handler)
        self.f_handler.setLevel(self.level)
        self.f_handler.setFormatter(self.f_format)


logger = Logger(str(Path(__file__).parent).split("/")[-1])