# © 2024 National Technology & Engineering Solutions of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.
# SPDX-License-Identifier: BSD-3-Clause

FROM python:3.9-slim

WORKDIR /app

COPY . .
RUN pip --proxy http://proxy.sandia.gov:80 --trusted-host pypi.org --trusted-host files.pythonhosted.org install -r requirements.txt

CMD ["python", "app.py"]
