# 
# Initial conditions for Conway's Game of Life: random, static, blinker and toad, pulsar, glieder, Gosper glider gun
# 
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
import matplotlib.animation as animation 

fig_size = (7, 7)
updateInterval = 50

#parameters
game_size = (70, 70)

##--------------------Initial Conditions------------------------------------
## random 
a = np.random.randint(0,2,size=game_size)
##-------------------------------------------------------------------------
## Static
#game_size = (6,21)
#a = np.zeros(game_size)
#a[2:4, 1:3] = 1
#a[1:4, 5:9] = [[0, 1, 1, 0],
#               [1, 0, 0, 1],
#               [0, 1, 1, 0]]
#a[1:5, 11:15] = [[0, 1, 1, 0],
#                 [1, 0, 0, 1],
#                 [0, 1, 0, 1],
#                 [0, 0, 1, 0]]
#a[1:4, 17:20] = [[1, 1, 0],
#                 [1, 0, 1],
#                 [0, 1, 0]]
##-------------------------------------------------------------------------
## Blinker
#game_size = (7,11)
#blinker = [1, 1, 1]
#toad = [[1, 1, 1, 0],
#        [0, 1, 1, 1]]
#a = np.zeros(game_size)
#a[2, 1:4] = blinker
#a[2:4, 6:10] = toad
##-------------------------------------------------------------------------
##Pulsar
#game_size = (17,17)
#a = np.zeros(game_size)
#a[2, 4:7] = 1
#a[4:7, 7] = 1
#a += a.T
#a += a[:, ::-1]
#a += a[::-1, :]
##-------------------------------------------------------------------------
## Glider
#game_size = (20,20)
#glider = [[1, 0, 0],
#          [0, 1, 1],
#          [1, 1, 0]]
#a = np.zeros(game_size)
#a[:3, :3] = glider
##-------------------------------------------------------------------------
# Gosper Glider Gun
#game_size = (70, 70)
#glider_gun =\
#[[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0],
 #[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0],
 #[0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1],
 #[0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1],
 #[1,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
 #[1,1,0,0,0,0,0,0,0,0,1,0,0,0,1,0,1,1,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0],
 #[0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0],
 #[0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],
 #[0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]]
#a = np.zeros(game_size)
#a[1:10,1:37] = glider_gun
##-------------------------------------------------------------------------

# filter: 3x3 matrix with ones except a 0 in the middle

def conway_mask(state, neighbor_num):
    # put in the state (0/1) of the relevant cell and the number of it's neighbors
    # returns whether cell will be reborn (1) or dies (or stays dead) (0)
    # number of alive neighbors:
    if(state==1):
        if(neighbor_num in [2,3]): return 1
        else: return 0
    elif(neighbor_num ==3): return 1
    else: return 0

def calc_grid(frame_num, img, table):
    # calculate one step of the game of life and make new image
    # frame_num is unused, but has to be passed for the animation library
    new_table = table.copy()
    filter = np.ones((3,3), dtype=int)
    filter[1][1] = 0
    neighbors = signal.convolve2d(table, filter, boundary="wrap", mode="same")
    for pos in [(x,y) for x in range(game_size[0]) for y in range(game_size[1])]:
        new_table[pos] = conway_mask(table[pos], neighbors[pos])
    
    # update data 
    img.set_data(new_table) 
    table[:] = new_table[:] 
    return img,

plt.imshow(a,interpolation='nearest',cmap=plt.cm.gray_r)
fig, ax = plt.subplots() 
img = ax.imshow(a, interpolation='nearest') 
ani = animation.FuncAnimation(fig, calc_grid, fargs=(img, a), 
                                  frames = 10, 
                                  interval=updateInterval, 
                                  save_count=50) 
plt.show()
# note: raises qt5ct errors depending on matplotlib version. can be ignored
