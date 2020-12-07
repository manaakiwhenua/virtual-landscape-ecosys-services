# Natcap InVEST Docker

This is published to https://hub.docker.com/r/richardlaw/natcap-invest, in lieu of an LCR Dockerhub account.

```bash
docker build -t richardlaw/natcap-invest:latest .
sudo docker push richardlaw/natcap-invest:latest
```

It follows instructions from https://invest.readthedocs.io/en/latest/installing.html, extrapolating them slightly since the instructions aren't complete. (I wish they maintained an official Docker build.)
