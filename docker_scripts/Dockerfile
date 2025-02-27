# syntax=docker/dockerfile:1

# Comments are provided throughout this file to help you get started.
# If you need more help, visit the Dockerfile reference guide at
# https://docs.docker.com/go/dockerfile-reference/

# Want to help us make this template better? Share your feedback here: https://forms.gle/ybq9Krt8jtBL3iCk7

FROM python:3.9

COPY requirements.txt .

RUN pip install --no-cache-dir -r requirements.txt
RUN apt-get update && apt-get install sudo
RUN sudo apt install ca-certificates

RUN wget -c --no-check-certificate https://mcell.org/downloads/Blender-2.93-CellBlender-4.0.6-Debian-GNU-Linux-9-20220708.tar.gz -O - | tar -xz
ENV CELLBLENDER_BASE_PATH="Blender-2.93-CellBlender/2.93/"
ENV MCELL_PATH="$CELLBLENDER_BASE_PATH/scripts/addons/cellblender/extensions/mcell/"
ENV PATH="$PATH:$CELLBLENDER_BASE_PATH/python/bin/"

COPY . . 

# Run the script command
CMD ["python","spatial_model.py"]
