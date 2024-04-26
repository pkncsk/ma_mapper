
import logging
#setup logger
logger = logging.getLogger('ma_mapper')
logger.setLevel(logging.DEBUG)
# File handler
file_handler = logging.FileHandler('ma_mapper.log')
file_handler.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s %(name)s %(levelname)s %(message)s')
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)

# Stream handler (for IDE/IPython/Jupyter)
stream_handler = logging.StreamHandler()
stream_handler.setLevel(logging.INFO)  # Adjust level as needed
stream_formatter = logging.Formatter('%(asctime)s %(levelname)s: %(message)s')
stream_handler.setFormatter(stream_formatter)
logger.addHandler(stream_handler)