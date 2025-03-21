import logging
import os
import datetime

# Define the logger globally
logger = None


def setup_logger():
    """
    Sets up a logging system where logs are stored in a separate folder and a single file per run.
    Ensures a single log file is used across multiple modules.
    """
    global logger
    if logger is not None:
        return logger  # Return existing logger to ensure single instance

    log_dir = os.path.join("workflows", "logs","logs")
    os.makedirs(log_dir, exist_ok=True)

    start_time = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    log_filename = f"Query_Logs_{start_time}.log"
    log_filepath = os.path.join(log_dir, log_filename)

    # Configure logging
    logger = logging.getLogger("linkedin_logger")
    logger.setLevel(logging.INFO)

    # Create file handler
    file_handler = logging.FileHandler(log_filepath)
    file_handler.setFormatter(logging.Formatter(
        '%(asctime)s - %(levelname)s - %(message)s'))

    # Create stream handler
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(logging.Formatter(
        '%(asctime)s - %(levelname)s - %(message)s'))

    # Add handlers to logger
    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)

    return logger
