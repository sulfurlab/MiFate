#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import xgboost as xgb

def read_table(path: str) -> pd.DataFrame:
    p = Path(path)
    if p.suffix.lower() in [".tsv", ".txt"]:
        return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    return pd.read_csv(path, dtype=str, keep_default_na=False)

def detect_id_col(df: pd.DataFrame) -> str:
    preferred = ["MAG", "Genome", "ID", "genome", "mag", "id", "sample", "Sample"]
    for c in preferred:
        if c in df.columns: return c
    return df.columns[0]

def normalize_cols(df: pd.DataFrame, id_col: str) -> pd.DataFrame:
    df = df.copy()
    new_cols = ["ID" if c == id_col else c.replace("KEGG_", "").upper() for c in df.columns]
    df.columns = new_cols
    return df

def load_features(path: str) -> list[str]:
    with open(path, "r", encoding="utf-8") as f:
        raw_feats = [x.strip().upper() for x in f if x.strip()]
        unique_feats = list(dict.fromkeys(raw_feats)) 
    return unique_feats

def load_thresholds(best_tsv: str):
    if not best_tsv or not Path(best_tsv).exists(): return 0.5, 0.5
    tab = pd.read_csv(best_tsv, sep='\t')
    sec = tab.loc[tab['model'].astype(str).str.upper() == 'SEC', 'threshold_youdenJ']
    dis = tab.loc[tab['model'].astype(str).str.upper() == 'DIS', 'threshold_youdenJ']
    return float(sec.values[0]) if len(sec) else 0.5, float(dis.values[0]) if len(dis) else 0.5

# -------------------------
# 核心：样本孤立化对齐
# -------------------------

def build_independent_matrix(raw_df: pd.DataFrame, target_features: list[str]):
    """
    不管输入文件有多少列，只提取模型需要的列。
    缺失的列直接视为该样本不持有该特征（补0）。
    """
    n_samples = len(raw_df)
    n_features = len(target_features)
    matrix = np.zeros((n_samples, n_features), dtype=np.float32)
    
    # 建立模型特征索引
    feat_to_idx = {feat: i for i, feat in enumerate(target_features)}
    
    # 仅遍历输入数据中存在的、且模型需要的列
    cols_to_process = [c for c in raw_df.columns if c in feat_to_idx]
    
    for col in cols_to_process:
        idx = feat_to_idx[col]
        # 转换为数值并二值化
        matrix[:, idx] = pd.to_numeric(raw_df[col], errors='coerce').fillna(0).gt(0).astype(np.float32)
            
    # 计算每个 MAG 命中模型特征的总数 (Individual Active Count)
    active_ones = matrix.sum(axis=1).astype(int)
    
    return matrix, active_ones

def classify_strictly_independent(probs, t, delta, active_counts, min_active):
    """
    分类逻辑不再使用全局 coverage。
    结果仅受：概率值、证据数量(active_counts)影响。
    """
    labels = np.full(len(probs), '', dtype=object)
    
    # 只有当这个 MAG 自身包含的功能基因太少时，才判定为 uncertain
    insufficient = (active_counts < min_active)
    # 概率在阈值附近的灰度带
    gray = (probs > t - delta) & (probs < t + delta)
    
    labels[insufficient | gray] = 'uncertain'
    labels[(labels == '') & (probs >= t + delta)] = 'hard'
    labels[(labels == '') & (probs <= t - delta)] = 'easy'
    
    # 剩余的根据阈值硬切
    mask = (labels == '')
    if mask.any():
        labels[mask] = np.where(probs[mask] >= t, 'hard', 'easy')
    return labels

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--infile', required=True)
    parser.add_argument('--sec_model', default='SEC_model.xgb')
    parser.add_argument('--dis_model', default='DIS_model.xgb')
    parser.add_argument('--sec_features', default='SEC_features.txt')
    parser.add_argument('--dis_features', default='DIS_features.txt')
    parser.add_argument('--best_thresholds', default='Best_thresholds.tsv')
    parser.add_argument('--out', default='MAG_pred_final.csv')
    parser.add_argument('--delta', type=float, default=0.05)
    parser.add_argument('--min_active', type=int, default=1, help="降低门槛防止全部变 uncertain")
    parser.add_argument('--chunk', type=int, default=5000)
    args = parser.parse_args()

    # 1. 加载
    sec_model = xgb.Booster(); sec_model.load_model(args.sec_model)
    dis_model = xgb.Booster(); dis_model.load_model(args.dis_model)
    sec_feats = load_features(args.sec_features)
    dis_feats = load_features(args.dis_features)
    t_sec, t_dis = load_thresholds(args.best_thresholds)

    # 2. 处理数据
    raw = read_table(args.infile)
    id_col = detect_id_col(raw)
    raw = normalize_cols(raw, id_col)
    ids = raw['ID'].values
    X_raw_only = raw.drop(columns=['ID'])

    # 3. 构建对齐矩阵 (每个样本的处理完全隔离)
    X_sec_np, active_sec = build_independent_matrix(X_raw_only, sec_feats)
    X_dis_np, active_dis = build_independent_matrix(X_raw_only, dis_feats)

    # 4. 预测
    dmat_sec = xgb.DMatrix(X_sec_np, feature_names=sec_feats)
    p_sec = sec_model.predict(dmat_sec)
    
    dmat_dis = xgb.DMatrix(X_dis_np, feature_names=dis_feats)
    p_dis = dis_model.predict(dmat_dis)

    # 5. 分类 (彻底去掉 min_coverage 参数)
    c_sec = classify_strictly_independent(p_sec, t_sec, args.delta, active_sec, args.min_active)
    c_dis = classify_strictly_independent(p_dis, t_dis, args.delta, active_dis, args.min_active)

    # 6. 导出
    pd.DataFrame({
        'ID': ids, 
        'pred_SEC': p_sec, 'pred_DIS': p_dis,
        'class_SEC': c_sec, 'class_DIS': c_dis,
        'active_genes_sec': active_sec,
        'active_genes_dis': active_dis
    }).to_csv(args.out, index=False)
    
    print(f"Success! Independent prediction saved to {args.out}")

if __name__ == '__main__':
    main()