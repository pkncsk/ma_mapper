print('test1')
import logging
#setup logger
logger = logging.getLogger('ma_mapper')
logger.setLevel(logging.DEBUG)
file_handler = logging.FileHandler('ma_mapper.log')
file_handler.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s %(name)s %(levelname)s %(message)s')
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)

logger.info('start logging, test')