# Madagascar Dockerfiles
You can pull images from [Docker Hub](https://hub.docker.com/r/madagascar/m8r).

For example:
```bash
$ docker pull madagascar/m8r:2.0-stable
```

## Running
```bash
docker run -it madagascar/m8r:2.0-stable /bin/bash
```

## Building
```bash
$ docker build -f ${RSFSRC}/admin/docker/[tag]/Dockerfile -t m8r .
```


