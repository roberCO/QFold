import numpy as np

def calculate_stats(data):

    mean_precisions = calculate_mean_precision(data)

    return mean_precisions

def calculate_mean_precision(data):

    minifold_tts = []
    minifold_precision = []

    random_tts = []
    random_precision = []

    for key in data.keys():

        if 'minifold' in key:

            minifold_tts.append(data[key]['min_tts'])
            minifold_precision.append(data[key]['precision'])

        if 'random' in key:

            random_tts.append(data[key]['min_tts'])
            random_precision.append(data[key]['precision'])

    mean_stats = {}
    mean_stats['minifold_precision'] = np.mean(minifold_precision)
    mean_stats['minifold_tts'] = np.mean(minifold_tts)
    mean_stats['random_precision'] = np.mean(random_precision)
    mean_stats['random_tts'] = np.mean(random_tts)

    return mean_stats
    
        