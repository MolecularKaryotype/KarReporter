import base64
from PIL import Image
import io


def reflect_image_horizontally_and_convert_to_base64(image_path, output_file):
    # Open an image file
    with Image.open(image_path) as img:
        # Reflect the image vertically
        reflected_image = img.transpose(Image.FLIP_LEFT_RIGHT)

        # Save the image to a bytes buffer
        buffered = io.BytesIO()
        reflected_image.save(buffered, format="PNG")

        # Encode the image to base64
        img_str = base64.b64encode(buffered.getvalue()).decode("utf-8")

    with open(output_file, 'w') as fp_write:
        fp_write.write(img_str)
    return img_str


# Example usage
image_path = 'magnifying_glass_icon.png'
base64_image = reflect_image_horizontally_and_convert_to_base64(image_path, 'magnifying_glass_icon_reflected.txt')
print(base64_image)