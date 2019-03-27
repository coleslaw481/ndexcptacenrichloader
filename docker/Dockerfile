# Use an official Python runtime as a parent image
FROM continuumio/miniconda3

RUN apt-get update
RUN apt-get install -y build-essential 

RUN mkdir /tmp/ndexcptacenrichloader
COPY ./ /tmp/ndexcptacenrichloader/
RUN pip install /tmp/ndexcptacenrichloader

RUN rm -rf /tmp/ndexcptacenrichloader

RUN mkdir /ndexcptacenrichloader
WORKDIR /ndexcptacenrichloader

CMD /bin/echo "To use run: ndexloadcptacenrich.py"
