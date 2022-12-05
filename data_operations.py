import numpy as np


def normalize_data(image: np.ndarray):
    image = np.nan_to_num(image)
    normalized_image = (image - image.min()) / (image.max() - image.min())

    return normalized_image


def make_sizes_match(image1, image2):
    # quit early if the image shapes already match
    if image1.shape == image2.shape:
        return image1, image2

    # image 1 is smaller
    if image1.size<image2.size:
        # resize image 1 to image2 (fill with zeros outter column/row)
        resized_1 = np.resize(image1, image2.shape)
        return resized_1, image2

    # image 2 is smaller
    elif image1.size>image2.size:
        # image 2 needs to be resized to image 1
        resized_2 = np.resize(image2, image1.shape)
        return image1, resized_2

    return None

