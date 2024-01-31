import os

from PIL import Image, ImageDraw, ImageFont


def crop_image(img, left, top, right, bottom):
    """
    Crop the image with the given margins.

    :param img: Image object to be cropped.
    :param left: Pixels to crop from the left.
    :param top: Pixels to crop from the top.
    :param right: Pixels to crop from the right.
    :param bottom: Pixels to crop from the bottom.
    :return: Cropped Image object.
    """
    width, height = img.size
    return img.crop((left, top, width - right, height - bottom))


def add_text_to_image(img, text):
    """
    Add text to the bottom-left corner of the image.

    :param img: Image object to add text to.
    :param text: Text to add.
    :return: Image object with text.
    """
    draw = ImageDraw.Draw(img)
    font = ImageFont.load_default()
    text_width, text_height = draw.textsize(text, font=font)

    position = (img.width - text_width - 10, img.height - text_height - 10)
    draw.text(position, text, font=font, fill="white")

    return img


def create_image_grid(images, grid_size):
    """
    Create a grid of images.

    :param images: List of Image objects.
    :param grid_size: Tuple (N, M) for grid dimensions.
    :return: Image object representing the grid.
    """
    img_width, img_height = images[0].size
    grid_img = Image.new("RGB", (img_width * grid_size[0], img_height * grid_size[1]))

    for y in range(grid_size[1]):
        for x in range(grid_size[0]):
            try:
                grid_img.paste(images[y * grid_size[0] + x], (x * img_width, y * img_height))
            except:
                pass

    return grid_img


# Example usage:
for condition in ['AF', 'SR']:
    image_folder = f"landmark_figures/{condition}"
    output_image_path = f"landmark_figures/landmark_{condition}.jpg"

    image_folder = f"ProbeViz/{condition}"
    output_image_path = f"ProbeViz/probes_{condition}.jpg"


    # Load all images from the folder
    image_files = [f for f in sorted(os.listdir(image_folder)) if os.path.isfile(os.path.join(image_folder, f))]
    images = [Image.open(os.path.join(image_folder, image_file)) for image_file in image_files]

    # Crop all images and add text
    cropped_and_texted_images = []
    dx = 200
    for idx, img in enumerate(images):
        cropped_img = crop_image(img, dx, dx, dx, dx)
        cropped_and_texted_images.append(cropped_img)

    N = 8
    M = 5
    # Create a grid of images (for example, 3x3)
    # NO CROPPING
    grid_image = create_image_grid(images, (N, M))
    # CROPING
    #grid_image = create_image_grid(cropped_and_texted_images, (N, M))
    grid_image.save(output_image_path)
