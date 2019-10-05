import sys
#ref  https://stackoverflow.com/questions/47733596/how-to-make-the-color-of-a-rectangle-change-without-doing-anything
import pygame

color_g = 0
color_r = 255

WHITE = (255, 255, 255)
size = (700, 500)

# --- start ---
pygame.init()
screen = pygame.display.set_mode(size)

# get current time
current_time = pygame.time.get_ticks()

# first change after 250 ms
change_color_time = current_time + 20

# --- mainloop ---
clock = pygame.time.Clock()

while True:

    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            sys.exit()

    # get current time 
    current_time = pygame.time.get_ticks()

    # check if it is time to change color
    if current_time >= change_color_time:
        # set new time to change color again
        change_color_time = current_time + 20
        # change color

        if color_r - 10  < 0:
            color_r = 0
            color_g += 5
            if color_g == 255:
                color_r = 255
                color_g = 0                    
        else:
            color_r -= 10

    screen.fill(WHITE)
    
    changing_Color = (color_r, color_g, 0)
    changing_C1 = (color_r, 97, 0)
    changing_C2 = (color_r, 165, 0)
    pygame.draw.rect(screen, changing_Color, [240,  60,  100,  40],0)
    
    pygame.draw.rect(screen, changing_C1,    [120,  60,  100,  40],0)
    pygame.draw.rect(screen, changing_C1,    [120, 200,  100,  40],0)
    
    pygame.draw.rect(screen, changing_C1,    [350, 200,  100,  40],0)
    pygame.draw.rect(screen, changing_C2,    [350,  60,  100, 100],0)
    pygame.draw.rect(screen, changing_C2,    [475,  60,  100,  40],0)
    pygame.display.flip()

    clock.tick(60)

pygame.quit()

