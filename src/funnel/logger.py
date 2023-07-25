import logging
import sys
import warnings

from loguru import logger

logger.remove()
logger.add(
    sys.stderr,
    format="|<blue>funnel</blue>|{level}| <green>{message}</green> ",
    colorize=True,
    level="DEBUG",
)

warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=UserWarning)
logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
