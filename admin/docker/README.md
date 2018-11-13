# Madagascar Dockerfiles
You can also pull images from [Docker Hub](https://hub.docker.com/r/madagascar/m8r).
For example,
```bash
$ docker pull madagascar/m8r:2.0-dev
```

## Building
```bash
$ docker build -f ./dockerfiles/[tag]/Dockerfile -t m8r .
```

## Running
```bash
docker run -it m8r /bin/bash
```
