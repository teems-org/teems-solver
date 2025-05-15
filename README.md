# TEEMS solver

[![License](https://img.shields.io/badge/License-GPL-blue.svg)](LICENSE)
[![Version](https://img.shields.io/badge/version-0.9-green.svg)](https://github.com/username/repo/releases)

This repository contains files necessary to build the optimization solver used within the R package TEEMS. Select HSL libraries are required and must be obtained directly from HSL (https://www.hsl.rl.ac.uk/).

## Table of Contents

- [Prerequisites](#prerequisites)
- [Installation](#installation)
  - [Docker installation](#docker-installation)
  - [HSL libraries](#hsl-libraries)
  - [Solver build](#solver-build)
    -[Layered build](#layered-build)
    -[Full build](#full-build)
  - [Cloning the Repository](#cloning-the-repository)
  - [Installing Dependencies](#installing-dependencies)
  - [Configuration](#configuration)
- [Usage](#usage)
- [Examples](#examples)
- [Contributing](#contributing)
- [Troubleshooting](#troubleshooting)
- [License](#license)
- [Contact](#contact)

## Prerequisites

The following must be installed/obtained prior to running the Docker build script:

- Docker https://www.docker.com/get-started/
- HSL library MA48 (version 2.2.0): https://www.hsl.rl.ac.uk/catalogue/ma48.html
- HSL library MA51 (version 1.0.0): https://www.hsl.rl.ac.uk/catalogue/ma51.html
- HSL library HSL_MC66 (version 2.2.0): https://www.hsl.rl.ac.uk/catalogue/hsl_mc66.html
- HSL library HSL_MP48 (version 2.1.1): https://www.hsl.rl.ac.uk/catalogue/hsl_mp48.html

## Installation

Instructions for installing the TEEMS solver.

### Docker installation
Follow the OS-specific installation instructions for the Docker containerization software: https://www.docker.com/get-started/

Linux users must ensure that Docker can be run without invoking sudo: https://docs.docker.com/engine/install/linux-postinstall/.

### HSL libraries
The required HSL libraries must be requested: https://www.hsl.rl.ac.uk
After recieving the tarballs for MA48, MA51, HSL_MC66, and HSL_MP48, copy the tarballs (e.g., ma48-2.2.0tar.gz) into the empty hsl folder at /teems-solver/hsl. Backware compatibility with previous HSL library versions is not guaranteed. This directory should now contain the following files:

- ma48-2.2.0.tar.gz
- ma51-1.0.0.tar.gz
- hsl_mc66-2.2.1.tar.gz
- hsl_mp48-2.1.1.tar.gz

### Solver build
In order to facilitate the solver build, a prebuilt Docker image with all dependencies is available.
The user must only link to their local HSL libraries to complete the build.
For developers and others preferring to build from scratch, a Dockerfile containing the full build is available as well.
The build time for the ``full build" is roughly 1 hour.

#### Prebuilt installation

#### Full build installation

Clone the repository to a local directory
Username is your git username (where you got the invite)
Password is: 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
```bash
git clone https://github.com/matthewcantele/teems-solver
```

Enter directory
```bash
cd teems-solver
```

Build the solver
Build time is approximately 1 hour depending on your local machine specs.
```bash

docker build --build-arg PATH_HSL_MA48="hsl/ma48-2.2.0.tar.gz" --build-arg PATH_HSL_MA51="hsl/ma51-1.0.0.tar.gz" --build-arg PATH_HSL_MC66="hsl/hsl_mc66-2.2.0.tar.gz" --build-arg PATH_HSL_MP48="hsl/hsl_mp48-2.1.1.tar.gz" -t teems:latest -f ./docker/full_build/Dockerfile .
```

Once built, check for the Docker image (teems:latest)
```bash
docker image ls
```

## Usage
The TEEMS solver is most easily utilized in conjunction with the TEEMS R package (link). It can however be called on solver-ready files. A middle ground option also exists with the in-situ-solve option within the TEEMS R package.

```bash
# Example command to run the software
./run.sh --option value

# Or for a library
import package_name

result = package_name.main_function()
```

## Examples

Include a few examples of how to use your software:

### Example 1: Basic Usage

```python
from package_name import feature

# Initialize
client = feature.Client()

# Use functionality
result = client.process_data("input")
print(result)
```

### Example 2: Advanced Configuration

```python
from package_name import feature

# Initialize with custom settings
config = {
    "option1": "value1",
    "option2": "value2"
}

client = feature.Client(config=config)
result = client.advanced_process("input", extra_param=True)
```

## Contributing

We welcome contributions to this project! Please follow these steps:

1. Fork the repository
2. Create a feature branch: `git checkout -b feature-name`
3. Commit your changes: `git commit -am 'Add some feature'`
4. Push to the branch: `git push origin feature-name`
5. Submit a pull request

Please read [CONTRIBUTING.md](CONTRIBUTING.md) for more details on our code of conduct and development process.

## Troubleshooting

### Common Issues

#### Issue 1: Error Message Description
Solution: Steps to resolve the issue.

#### Issue 2: Error Message Description
Solution: Steps to resolve the issue.

For more help, please check our [FAQ](docs/FAQ.md) or [open an issue](https://github.com/username/repository/issues).

## License

This project is licensed under the GPLv3.0 License - see the [LICENSE](LICENSE) file for details.

## Code authorship
This work is the culmination of many years of efforts and collaborations. The C source code (src) and main build was written by Tom Kompas and Ha Van Pham. The binary parsing code (bin_parser) was contributed by Martin Ingrahm. Finally, Matthew Cantele authored the Docker and Singularity scripts.

## Contact

- Project Maintainer: [Matthew Cantele](mailto:matthew.cantele@protonmail.com)
- Project Homepage: [https://github.com/matthewcantele/teems-solver](https://github.com/matthewcantele/teems-solver)
- Bug Reports: [https://github.com/username/repository/issues](https://github.com/username/repository/issues)