#! /usr/bin/env bash

isort src tests

black src tests

flake8 src tests