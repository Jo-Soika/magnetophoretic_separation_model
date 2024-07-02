def hit_bottom(time, particle_position):
    return particle_position[2]


def travel_too_far(time, particle_position, max_distance):
    return particle_position[0] - max_distance


def hit_upper_wall(time, particle_position, channel_height):
    return particle_position[2] - channel_height
