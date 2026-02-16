from typing import Literal
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from rich import print


def overlap_poincare(
    folder,
    mode: Literal["_reflection", "_transmission", ""],
    version: str,
    s1_col=" Normalized s 1 ",
    s2_col=" Normalized s 2 ",
    s3_col=" Normalized s 3 ",
):

    df_h_mode = pd.read_csv(f"./{folder}/H{mode}{version}.csv", sep=";", skiprows=8)
    df_v_mode = pd.read_csv(f"./{folder}/V{mode}{version}.csv", sep=";", skiprows=8)
    df_a_mode = pd.read_csv(f"./{folder}/A{mode}{version}.csv", sep=";", skiprows=8)
    df_d_mode = pd.read_csv(f"./{folder}/D{mode}{version}.csv", sep=";", skiprows=8)
    df_r_mode = pd.read_csv(f"./{folder}/R{mode}{version}.csv", sep=";", skiprows=8)
    df_l_mode = pd.read_csv(f"./{folder}/L{mode}{version}.csv", sep=";", skiprows=8)

    S1_h = df_h_mode[s1_col].mean()
    S2_h = df_h_mode[s2_col].mean()
    S3_h = df_h_mode[s3_col].mean()
    s_h = [S1_h, S2_h, S3_h]

    S1_v = df_v_mode[s1_col].mean()
    S2_v = df_v_mode[s2_col].mean()
    S3_v = df_v_mode[s3_col].mean()
    s_v = [S1_v, S2_v, S3_v]

    S1_a = df_a_mode[s1_col].mean()
    S2_a = df_a_mode[s2_col].mean()
    S3_a = df_a_mode[s3_col].mean()
    s_a = [S1_a, S2_a, S3_a]

    S1_d = df_d_mode[s1_col].mean()
    S2_d = df_d_mode[s2_col].mean()
    S3_d = df_d_mode[s3_col].mean()
    s_d = [S1_d, S2_d, S3_d]

    S1_r = df_r_mode[s1_col].mean()
    S2_r = df_r_mode[s2_col].mean()
    S3_r = df_r_mode[s3_col].mean()
    s_r = [S1_r, S2_r, S3_r]

    S1_l = df_l_mode[s1_col].mean()
    S2_l = df_l_mode[s2_col].mean()
    S3_l = df_l_mode[s3_col].mean()
    s_l = [S1_l, S2_l, S3_l]

    s_h_ideal = [1, 0, 0]
    s_v_ideal = [-1, 0, 0]
    s_a_ideal = [0, -1, 0]
    s_d_ideal = [0, 1, 0]
    s_r_ideal = [0, 0, 1]
    s_l_ideal = [0, 0, -1]

    overlap_h = overlap_polarizations(s_h, s_h_ideal)
    overlap_v = overlap_polarizations(s_v, s_v_ideal)
    overlap_a = overlap_polarizations(s_a, s_d_ideal)
    overlap_d = overlap_polarizations(s_d, s_a_ideal)
    overlap_r = overlap_polarizations(s_r, s_l_ideal)
    overlap_l = overlap_polarizations(s_l, s_r_ideal)

    avg_overlap = (
        overlap_h + overlap_v + overlap_a + overlap_d + overlap_r + overlap_l
    ) / 6

    print(f"Overlap H: {overlap_h}")
    print(f"Overlap V: {overlap_v}")
    print(f"Overlap A: {overlap_a}")
    print(f"Overlap D: {overlap_d}")
    print(f"Overlap R: {overlap_r}")
    print(f"Overlap L: {overlap_l}")

    print("-----------------------")
    print(f"Average Overlap: {avg_overlap}")
    print("-----------------------")


def overlap_polarizations(s, p):
    return 1 / 2 * (1 + np.dot(s, p))


overlap_poincare(folder="./data/", mode="_reflection", version="")
