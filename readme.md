# Docker crypto

Build
> docker build https://raw.githubusercontent.com/Mister7F/docker-crypto/master/debian/Dockerfile -t "crypto-debian"

> docker build https://raw.githubusercontent.com/Mister7F/docker-crypto/master/sage_10/Dockerfile -t "crypto-sage_10"

Start
> docker run -it --rm crypto-debian:latest

> docker run -it --rm crypto-sage_10:latest

All the code are not from me, I just added scripts to be able to call them in the CLI.

# Sources
- https://github.com/WardBeullens/BreakingRainbow (faster with lower sage version for some reason)
- https://github.com/GiacomoPope/Castryck-Decru-SageMath
