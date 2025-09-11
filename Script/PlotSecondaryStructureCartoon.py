#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 19:55:44 2025

@author: mayank
"""
import sys
# conda_env_path = '/var/www/miniconda3/envs/biotite/lib/python3.13/site-packages'
# sys.path.insert(0, conda_env_path)
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle
import biotite
import biotite.sequence as seq
import biotite.sequence.graphics as graphics

# Create 'FeaturePlotter' subclasses
# for drawing the scondary structure features


class HelixPlotter(graphics.FeaturePlotter):
    def __init__(self):
        pass

    # Check whether this class is applicable for drawing a feature
    def matches(self, feature):
        if feature.key == "SecStr":
            if "sec_str_type" in feature.qual:
                if feature.qual["sec_str_type"] == "helix":
                    return True
        return False

    # The drawing function itself
    def draw(self, axes, feature, bbox, loc, style_param):
        # Approx. 1 turn per 3.6 residues to resemble natural helix
        n_turns = np.ceil((loc.last - loc.first + 1) / 3.6)
        x_val = np.linspace(0, n_turns * 2 * np.pi, 100)
        # Curve ranges from 0.3 to 0.7
        y_val = (-0.4 * np.sin(x_val) + 1) / 2

        # Transform values for correct location in feature map
        x_val *= bbox.width / (n_turns * 2 * np.pi)
        x_val += bbox.x0
        y_val *= bbox.height
        y_val += bbox.y0

        # Draw white background to overlay the guiding line
        background = Rectangle(
            bbox.p0, bbox.width, bbox.height, color="white", linewidth=0
        )
        axes.add_patch(background)
        axes.plot(x_val, y_val, linewidth=2, color='#ff4d6d')
                  # color=biotite.colors["dimgreen"])


class SheetPlotter(graphics.FeaturePlotter):
    def __init__(self, head_width=0.8, tail_width=0.5):
        self._head_width = head_width
        self._tail_width = tail_width

    def matches(self, feature):
        if feature.key == "SecStr":
            if "sec_str_type" in feature.qual:
                if feature.qual["sec_str_type"] == "sheet":
                    return True
        return False

    def draw(self, axes, feature, bbox, loc, style_param):
        x = bbox.x0
        y = bbox.y0 + bbox.height / 2
        dx = bbox.width
        dy = 0

        if loc.defect & seq.Location.Defect.MISS_RIGHT:
            # If the feature extends into the prevoius or next line
            # do not draw an arrow head
            draw_head = False
        else:
            draw_head = True

        axes.add_patch(
            biotite.AdaptiveFancyArrow(
                x,
                y,
                dx,
                dy,
                self._tail_width * bbox.height,
                self._head_width * bbox.height,
                # Create head with 90 degrees tip
                # -> head width/length ratio = 1/2
                head_ratio=0.5,
                draw_head=draw_head,
                # color=biotite.colors["orange"],
                color='#ffc600',
                linewidth=0,
            )
        )
# 🔹 Convert Secondary Structure String to Biotite Annotation
def sec_str_to_annotation(sec_str):
    features = []
    i = 0
    while i < len(sec_str):
        start = i
        if sec_str[i] == "H":  # Helix
            while i < len(sec_str) and sec_str[i] == "H":
                i += 1
            features.append(seq.Feature("SecStr", [seq.Location(start, i)], {"sec_str_type": "helix"}))
        elif sec_str[i] == "E":  # Sheet
            while i < len(sec_str) and sec_str[i] == "E":
                i += 1
            features.append(seq.Feature("SecStr", [seq.Location(start, i)], {"sec_str_type": "sheet"}))
        i += 1
    return seq.Annotation(features)

def PlotSecondaryStructureCartoon(sec_str,JobID,ID):
    # 🟢 **Secondary Structure Input String**
    
    
    # Convert to Biotite annotation
    annotation = sec_str_to_annotation(sec_str)
    
    # **Plot**
    fig = plt.figure(figsize=(8.0, 3.0))
    ax = fig.add_subplot(111)
    graphics.plot_feature_map(
        ax,
        annotation,
        symbols_per_line=60,
        show_numbers=True,
        show_line_position=True,
        loc_range=(1, len(sec_str) + 1),
        feature_plotters=[HelixPlotter(), SheetPlotter()],
    )
    fig.tight_layout()
    plt.savefig("Jobs/"+JobID+"/"+str(ID)+".png", dpi=300)
# sec_str = "LLLLLHHHHHHLLLLLEEEEELLLLLHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHLLLLLEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEELLLLLHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL"
# JobID='test'
# PlotSecondaryStructureCartoon(sec_str,JobID)
if __name__ == '__main__':
    globals()[sys.argv[1]](sys.argv[2],sys.argv[3],sys.argv[4])