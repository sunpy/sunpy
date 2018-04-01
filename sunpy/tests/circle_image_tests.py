"""
A script to check the image test results on CircleCI, and upload
any new images to the sunpy-figure-tests repository.
"""
import os

image_directory = None
for f in os.listdir('.'):
    if 'result_images' in f:
        image_directory = f
        print(image_directory)

if image_directory is None:
    raise RuntimeError("Could not find directory with saved images")
