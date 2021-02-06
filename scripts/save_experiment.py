#!/usr/bin/env python
# coding: utf-8

import os
import shutil
from datetime import date

exp_files = ["evolution_example.txt", "configs/sample.ini"]

MAIN_OUTPUT = "save_initial_tests"


def save_files(name):
    today = date.today()
    d1 = today.strftime("%d%m%Y")
    name = d1 + "_" + name
    name = MAIN_OUTPUT + "/" + name
    os.makedirs(name, exist_ok=True)
    for f in exp_files:
        fname = f if "/" not in f else f.split("/")[-1]
        shutil.copyfile(f, name + "/" + fname)


def main():
    name = input("name > ")
    if len(name) == "0":
        print("error name")
        quit()
    save_files(name)


if __name__ == "__main__":
    main()
