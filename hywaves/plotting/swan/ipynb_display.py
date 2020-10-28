#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd
import ipywidgets as widgets
from IPython import display


def show_table(df):
    # create output widgets
    widget = widgets.Output()
    # render in output widgets
    with widget:
        with pd.option_context("display.max_rows", None, "display.max_columns", None):
            display.display(df)


