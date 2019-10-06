import pygame
import sys
import random
import numpy as np
from pygame.locals import *

pygame.init()

screen = pygame.display.set_mode((800, 640))
#ref http://adamdempsey90.github.io/python/bouncing_ball/bouncing_ball.html
#  https://www.101computing.net/pong-tutorial-using-pygame-adding-a-bouncing-ball/
pygame.display.set_caption("Simulate bouncing ball")
gravity = 4.5
vel1 = 0
vel2 = 0
Y1 = random.randint(20, 50)
Y2 = random.randint(40,100)

all_s1 = []
all_s2 = []

while True:
    for event in pygame.event.get():
        if event.type in (QUIT, KEYDOWN):
            sys.exit()
    screen.fill((255, 255, 255))

    vel1 += gravity
   
    if Y1 >= 610:
        vel1 = int( -vel1 *0.85 ) 
    Y1  +=  int(vel1)
    
    vel2 += gravity

    if Y2 >= 610:
        vel2 = int( -vel2 *0.80 ) 
    Y2  +=  int(vel2)

    all_s1.append(Y1)
    all_s2.append(Y2)

    if len(all_s1) >= 200 and abs(np.mean(all_s1[100:]) - np.median(all_s1[100:])) < 2.5:
        all_s1 = all_s1[-200:]
        Y1 = 610
    if len(all_s2) >= 200 and abs(np.mean(all_s2[100:]) - np.median(all_s2[100:])) < 2:
        all_s2 = all_s2[-200:]
        Y2 = 610
    pygame.draw.circle(screen, (220, 250, 30 ), (200, Y1) , 30, 0)
    pygame.draw.circle(screen, (233, 150, 122), (450, Y2) , 30, 0)

    pygame.display.update()

