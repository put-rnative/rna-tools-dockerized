FROM python:3.10

# This commit id is known to work
COPY lociPARSE-8c7acbe4e7c486122a4c261b1ea68fff7247b796 /opt/lociPARSE

RUN cd /opt/lociPARSE \
 && pip install .

CMD ["/bin/bash"]

COPY lociPARSE /usr/bin/lociPARSE
