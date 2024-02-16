from ctypes.wintypes import tagRECT
from ssl import CHANNEL_BINDING_TYPES
import pygame
import random
import time
from dataclasses import dataclass
#import pylsl
from pylsl import StreamInfo, StreamOutlet

<<<<<<< HEAD
# pylsl config
# Define the StreamInfo object
info = StreamInfo(name='game', channel_count=2,channel_format='string')
# Create the StreamOutlet object
outlet = StreamOutlet(info, chunk_size=0, max_buffered=360)
# Destroy the StreamInfo object to save space (optional)
info.__del__()
# Check whether there are consumers connected to the outlet
outlet.have_consumers()
# Wait for consumers to show up without wasting resources
# Returns True if the wait was successful, False if the timeout expired
# Only turns true if inlet.open_stream() was called
outlet.wait_for_consumers(timeout=10)

# pygame

=======
>>>>>>> accee9b5c61a779940a3b83427da040bb17d2f99
clock = pygame.time.Clock()


class Game:
    def __init__(self) -> None:
        self.number_of_experiments: int = 0
        self.white: any = (255, 255, 255)
        self.black: any  = (0, 0, 0)
        self.grey: any = (200, 200, 200)
        self.red: any  = (255, 0, 0)
        self.blue: any = (0, 0, 255)
        self.green: any = (0, 255, 0)
        self.box_reached: bool  = False
        self.chosen_box: any = []
        # Define the marker stream
        info = StreamInfo(name='pygame_markers',
                          type='Markers',
                          channel_count=1,
                          nominal_srate=0,
                          channel_format='string',
                          source_id='pygame_markers')
        self.marker_stream = StreamOutlet(info)
        pygame.init()
        self.options = [[175, 275, 50, 50], [575, 275, 50, 50], [375, 75, 50, 50]]
        self.dis = pygame.display.set_mode((800, 500))
        pygame.display.set_caption('BCI Experiment')
        self.initialize_canvas()
        self.keyboard_controls()

    def initialize_canvas(self):
        self.dis.fill(self.white)
        if self.chosen_box:
            draw_box = [self.chosen_box[0]-10, self.chosen_box[1]-10, self.chosen_box[2]+20, self.chosen_box[3]+20]
            pygame.draw.rect(self.dis, self.grey, draw_box)
        for box in self.options:
            pygame.draw.rect(self.dis, self.blue, box)

        pygame.display.update()

    def reset_game(self, timeout=False):
        self.number_of_experiments += 1
        self.chosen_box = random.choice(self.options)
        self.dis.fill(self.white)
        pygame.display.update()
        if timeout:
            self.display_break()
        time.sleep(1)
        self.box_reached = False
        pygame.draw.rect(self.dis, self.blue, [375, 275, 50, 50])

    def display_break(self):
        keyboard_input = False
        self.dis.fill(self.white)
        font = pygame.font.Font('freesansbold.ttf', 32)
        text = font.render('Press any Key to continue', True, self.grey, self.white)
        textRect = text.get_rect()
        self.dis.blit(text, textRect)
        pygame.display.update()
        while not keyboard_input:
            for event in pygame.event.get():
                if event.type == pygame.KEYDOWN:
                    keyboard_input = True


    def keyboard_controls(self):
        self.reset_game()
        if self.number_of_experiments % 10 == 0:
            print(self.number_of_experiments)
            print("Going to a break")
            self.display_break()

        self.initialize_canvas()
        
        target_set = False

        while not target_set:
            for event in pygame.event.get():
                    if event.type == pygame.QUIT:
                        target_set = True
                    if event.type == pygame.KEYDOWN:
                        if event.key == pygame.K_LEFT:
                            target_set = True
                            self.move_box((-5, 0))
                        elif event.key == pygame.K_RIGHT:
                            target_set = True
                            self.move_box((5, 0))
                        elif event.key == pygame.K_UP:
                            target_set = True
                            self.move_box((0, -5))

    def move_box(self, direction):
        if random.uniform(0, 1) < 0.3:
            print("richtung messed up")
            possible_directions = [(-5, 0), (5, 0), (0, -5)]
            possible_directions.remove(direction)
            direction = random.choice(possible_directions)

        x1 = 375
        y1 = 275

        while not self.box_reached:
            x1 += direction[0]
            y1 += direction[1]
            self.initialize_canvas()
            pygame.draw.rect(self.dis, self.black, [x1, y1, 50, 50])
            pygame.display.update()
            if x1 == 175 and y1 == 275:
                self.box_reached = True
                self.change_color(x1, y1)
            elif x1 == 575 and y1 == 275:
                self.box_reached = True
                self.change_color(x1, y1)
            elif x1 == 375 and y1 == 75:
                self.box_reached = True
                self.change_color(x1, y1)
            
            clock.tick(60)

    def change_color(self, x, y):
        if x ==  self.chosen_box[0] and y == self.chosen_box[1]:
            color = self.green
<<<<<<< HEAD
            outlet.push_sample("S1")
=======
            self.marker_stream.push_sample(['S1'])
>>>>>>> accee9b5c61a779940a3b83427da040bb17d2f99
            pygame.mixer.music.load('ding.mp3')
            pygame.mixer.music.play(0)
        else:
            color = self.red
<<<<<<< HEAD
            outlet.push_sample("S2")
=======
            self.marker_stream.push_sample(['S2'])
>>>>>>> accee9b5c61a779940a3b83427da040bb17d2f99
        pygame.draw.rect(self.dis, color, [x-10, y-10, 70, 70])
        pygame.draw.rect(self.dis, self.blue, [x, y, 50, 50])
        pygame.display.update()
        time.sleep(1)
        self.keyboard_controls()

lets_play = Game()

