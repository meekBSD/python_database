import pygame
from pygame.locals import *
from sys import exit

import math

pygame.init()
screen = pygame.display.set_mode((640,500))
pygame.display.set_caption("Heart Line - pygame window") 

# We need this if we want to be able to specify our
#  arc in degrees instead of radians

def degreesToRadians(deg):
    return deg/180.0 * math.pi


# Draw an arc that is a portion of a circle.
# We pass in screen and color,
# followed by a tuple (x,y) that is the center of the circle, and the radius.
# Next comes the start and ending angle on the "unit circle" (0 to 360)
#  of the circle we want to draw, and finally the thickness in pixels

def drawCircleArc(screen,color,center,radius,startDeg,endDeg,thickness):
    (x,y) = center
    rect = (x-radius,y-radius,radius*2,radius*2)
    startRad = degreesToRadians(startDeg)
    endRad = degreesToRadians(endDeg)
   
    pygame.draw.arc(screen,color,rect,startRad,endRad,thickness)
    

white = (255,255,255);
red = (255,0,0);
green = (0,255,0);
blue = (0,0,255);

while True:

    for event in pygame.event.get():
        if event.type == QUIT:
            pygame.quit(); exit();


    background = pygame.Surface(screen.get_size())
    background = background.convert()
    background.fill((255, 255, 255))
    screen.blit(background, (0,0))

    # Part of a red circle, arc from 0 to 180 degrees
    # Center is at 200,150, and radius is 120
    drawCircleArc(screen,red,(200,150),120,0,180,3)

    # Part of a blue circle, arc from 0 to 180 degrees
    # Center is at 440,150, radius is 120, thickness is 2
    drawCircleArc(screen, red,(440,150), 120, 0,180,3)

    
    # Part of a green circle, arc from -48 to 3 degrees
    # Center is at 320,150, radius is 240, thickness is 3
    drawCircleArc(screen,red,(320,150),240,-46,3,3)
    drawCircleArc(screen,red,(320,150),240,177,230,3)
    #drawCircleArc(screen,blue,(300,150),100,-45,+45,3)

    #pygame.draw.aaline(screen, (0, 0, 0), (480, 425), (550, 325), 1)
    pygame.draw.line(screen, red, (320 - 240/math.sqrt(2), 150+240/math.sqrt(2)), (320, 150 + 240 * math.sqrt(2)), 3)
    pygame.draw.line(screen, red, (320+240/math.sqrt(2), 150+240/math.sqrt(2)), (320, 150 + 240 * math.sqrt(2)), 3)

    #screen.blit(background, (0,0))
    pygame.display.flip()
    #pygame.display.update()

