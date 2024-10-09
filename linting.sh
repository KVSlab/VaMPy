#! /usr/bin/env bash

flake8 src tests

isort src tests

black src tests

mypy src tests
