# dataset
TRAIN_DATASET = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '12', '13', '14']
TEST_DATASET = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '12', '13', '14']
FRAMES = 2
VIEWS = 4

# training
TRAIN_GPUS = '0'
EX_ID = 'model'
LOSS = ['pairwise', 'triplewise']
LEARNING_RATE = 1e-5
MAX_EPOCH = 1000

NETWORK = 'resnet'
RE_ID = '0'
if RE_ID:
    TRAIN_RESUME = "./models/" + RE_ID + '.pth'
MODEL_SAVE_NAME = "./models/" + EX_ID + '.pth'

# parameters
MARGIN = 0.5
DATASET_SHUFFLE = 0
LOADER_SHUFFLE = 1

# inference
INF_ID = 'model'
DISPLAY = 0
INFTY_COST = 1e+5
RENEW_TIME = 30

