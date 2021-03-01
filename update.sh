#!/bin/bash

mongoexport -d TOV-LE -c Exponents-polytropic-EOS -o Aspects/database/database.txt
pipreqs .
