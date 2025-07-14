#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MetabarcodingTaxonomy 클래스:
- ASV 필터링 및 잘라내기
- Unclassified 비율 계산 및 시각화
- Retained taxa count 비율 시각화
- QIIME2 스타일 누적(스택) 바플롯 (level-*_truncated.csv → level-N_barplot.pdf)

출력 파일은 각 입력 파일이 있던 동일 디렉토리에 PDF로 저장됩니다.
"""
import re
import os
import glob
import datetime
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

# 전체 폰트 Arial
plt.rcParams['font.family'] = 'Arial'

class MetabarcodingTaxonomy:
    DEFAULT_LEVEL_LABELS = [
        'Kingdom', 'Phylum', 'Class',
        'Order', 'Family', 'Genus', 'Species'
    ]

    def __init__(
        self,
        input_dir='.',
        level_pattern='level-*.csv',
        sample_col=None,
        level_labels=None
    ):
        self.input_dir = input_dir
        pattern = os.path.join(input_dir, level_pattern)
        all_paths = sorted(glob.glob(pattern))
        # 원본 level-N.csv만 필터링
        self.file_paths = [
            fp for fp in all_paths
            if re.match(r'^level-\d+\.csv$', os.path.basename(fp))
        ]
        if not self.file_paths:
            raise FileNotFoundError(f"No level files found in {input_dir}")
        self.level_names = [
            os.path.splitext(os.path.basename(fp))[0]
            for fp in self.file_paths
        ]
        # X축 레이블 세팅
        if level_labels:
            if len(level_labels) != len(self.level_names):
                raise ValueError("level_labels length must match number of levels")
            self.level_labels = level_labels
        else:
            self.level_labels = MetabarcodingTaxonomy.DEFAULT_LEVEL_LABELS
        self.sample_col = sample_col
        self.stats_df = None
        self.taxa_count_df = None
        self.log(f"Initialized with levels: {self.level_names}")

    def log(self, message: str):
        now = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        print(f"[{now}] {message}")

    def is_filtered_taxon(self, taxon: str) -> bool:
        segs = taxon.split(';')
        lower = taxon.lower()
        # Incertae 포함 또는 _sp
        if 'incertae' in lower or taxon.endswith('_sp'):
            return True
        # phylum 이후 모두 '__'
        if len(segs) >= 3 and all(s == '__' for s in segs[2:]):
            return True
        # kingdom 이후 모두 '__'
        if len(segs) >= 2 and all(s == '__' for s in segs[1:]):
            return True
        # 숫자나 '-' 포함
        if re.search(r'[0-9\-]', taxon):
            return True
        # 기타 키워드 필터링
        if any(k in lower for k in (
            'uncultured','unidentified','candidum',
            'candidatus','metagenome'
        )):
            return True
        # 특수문자
        if re.search(r"[^A-Za-z0-9_;]", taxon):
            return True
        return False

    @staticmethod
    def _last_tax_label(tax_str: str) -> str:
        if tax_str == 'Other':
            return 'Other'
        parts = tax_str.split(';')
        last = parts[-1]
        return last.split('__',1)[1] if '__' in last else last

    def filter_and_truncate(self, df: pd.DataFrame, source_path: str):
        self.log(f"Processing {source_path}")
        sample_col = self.sample_col or df.columns[0]
        taxa_cols = df.columns.drop(sample_col)
        filtered = [c for c in taxa_cols if self.is_filtered_taxon(c)]
        retained = [c for c in taxa_cols if c not in filtered]
        self.log(f" Columns total={len(taxa_cols)}, filtered={len(filtered)}, retained={len(retained)}")

        def truncate(tok):
            parts = tok.split(';')
            kept = []
            for p in parts:
                lp = p.lower()
                if (
                    p=='__' or 'incertae' in lp or p.endswith('_sp') or
                    re.search(r'[0-9\-]', p) or
                    any(k in lp for k in (
                        'uncultured','unidentified','candidum',
                        'candidatus','metagenome'
                    )) or
                    re.search(r"[^A-Za-z0-9_;]", p)
                ):
                    break
                kept.append(p)
            return ';'.join(kept)

        df_filtered = df[[sample_col] + filtered]
        df_retained = df[[sample_col] + retained]
        df_trunc = df.copy()
        df_trunc.columns = [sample_col] + [truncate(c) for c in taxa_cols]

        total = df[taxa_cols].sum().sum()
        fsum  = df_filtered[filtered].sum().sum() if filtered else 0
        pct   = fsum / total * 100 if total > 0 else 0
        self.log(f" Total count={total}, filtered count={fsum} ({pct:.2f}%)")

        base = os.path.splitext(os.path.basename(source_path))[0]
        dirp = os.path.dirname(source_path)
        df_filtered.to_csv(f"{dirp}/{base}_filtered.csv", index=False)
        df_retained.to_csv(f"{dirp}/{base}_retained.csv", index=False)
        df_trunc.to_csv(f"{dirp}/{base}_truncated.csv", index=False)

    def compute_unclassified_stats(self):
        self.log("Computing unclassified stats...")
        stats = {}
        samples = None
        for fp, lvl in zip(self.file_paths, self.level_names):
            df = pd.read_csv(fp)
            sc = self.sample_col or df.columns[0]
            if samples is None:
                samples = df[sc].tolist()
            taxa = df.columns.drop(sc)
            total = df[taxa].sum(axis=1)
            filt_cols = [c for c in taxa if self.is_filtered_taxon(c)]
            fsum = df[filt_cols].sum(axis=1) if filt_cols else pd.Series(0, index=df.index)
            stats[lvl] = (fsum / total.replace(0, pd.NA)).values
            self.log(f" Level {lvl}: {len(df)} samples")
        self.stats_df = pd.DataFrame(stats, index=samples)
        return self.stats_df

    def compute_taxa_counts(self):
        self.log("Computing retained taxa count ratios...")
        ratios = {}
        for fp, lvl in zip(self.file_paths, self.level_names):
            df = pd.read_csv(fp)
            sc = self.sample_col or df.columns[0]
            taxa = df.columns.drop(sc)
            total_cols = len(taxa)
            retained = len([c for c in taxa if not self.is_filtered_taxon(c)])
            ratios[lvl] = retained / total_cols
            self.log(f" Level {lvl}: retained {retained}/{total_cols}")
        self.taxa_count_df = pd.Series(ratios, name='retained_taxa_ratio')
        return self.taxa_count_df

    def plot_well_classified(self):
        self.log("Plotting well-classified proportions...")
        df = self.stats_df
        well = 1 - df
        fig, ax = plt.subplots(figsize=(9,5), dpi=450)
        for samp in well.index:
            ax.plot(
                self.level_labels,
                well.loc[samp] * 100,
                marker='o', linewidth=4, alpha=0.8, label=samp, markersize=10
            )
        ax.set_ylim(0, 105)
        ax.grid(True, linestyle='--', alpha=0.4)
        ax.legend(loc='best', fontsize=20)
        ax.set_ylabel('Well-classified (%)', fontsize=20)
        ax.set_xticklabels(self.level_labels, fontsize=20)
        # Y축 tick label size
        ax.tick_params(axis='y', labelsize=15)
        plt.tight_layout()
        out = os.path.join(self.input_dir, 'Well_classified_proportion.pdf')
        fig.savefig(out, dpi=450, format='pdf')
        self.log(f"Saved {out}")
        plt.close(fig)

    def plot_taxa_retained(self):
        self.log("Plotting retained taxa counts...")
        df = self.compute_taxa_counts()
        fig, ax = plt.subplots(figsize=(9,5), dpi=450)
        ax.plot(
            self.level_labels,
            df.values * 100,
            marker='o', linewidth=4, markersize=10
        )
        ax.set_ylim(0, 105)
        ax.grid(True, linestyle='--', alpha=0.4)
        ax.set_ylabel('Retained taxa (%)', fontsize=20)
        ax.set_xticklabels(self.level_labels, fontsize=20)
        plt.tight_layout()
        out = os.path.join(self.input_dir, 'Retained_taxa_ratio.pdf')
        fig.savefig(out, dpi=450, format='pdf')
        self.log(f"Saved {out}")
        plt.close(fig)

    def plot_cumulative_barplots(self, dpi: int = 300, top_n: int = 10):
        self.log("Plotting cumulative barplots for each level...")
        pattern = os.path.join(self.input_dir, 'level-*_truncated.csv')
        t_files = sorted(glob.glob(pattern))
        if not t_files:
            raise FileNotFoundError("No truncated files to plot.")

        for fp in t_files:
            level = os.path.basename(fp).split('_')[0]
            idx = int(level.split('-')[1]) - 1
            rank = (
                self.level_labels[idx]
                if 0 <= idx < len(self.level_labels)
                else level
            )

            df = pd.read_csv(fp)
            sc = df.columns[0]
            counts = df.set_index(sc)
            rel = counts.div(counts.sum(axis=1), axis=0) * 100

            means = rel.mean(axis=0).sort_values(ascending=False)
            top_labels = means.index[:top_n].tolist()

            top_rel = rel[top_labels].copy()
            other_val = (100 - top_rel.sum(axis=1)).clip(lower=0)

            rel2 = top_rel.copy()
            rel2['Other'] = other_val
            cols = top_labels + ['Other']

            cmap = plt.get_cmap('Set1')
            colors = [cmap(i) for i in range(top_n)] + ['black']

            fig, ax = plt.subplots(figsize=(8,5), dpi=dpi)
            bottom = np.zeros(len(rel2))
            x = np.arange(len(rel2))
            for col, color in zip(cols, colors):
                ax.bar(
                    x, rel2[col], bottom=bottom,
                    label=MetabarcodingTaxonomy._last_tax_label(col),
                    color=color, edgecolor='none'
                )
                bottom += rel2[col].values

            ax.set_xticks(x)
            ax.set_xticklabels(rel2.index, rotation=30, ha='right', fontsize=8)
            ax.set_ylim(0, 100)
            ax.set_ylabel('Relative abundance (%)', fontsize=20)
            ax.grid(axis='y', linestyle='--', alpha=0.3)

            ax.legend(
                loc='upper left',
                bbox_to_anchor=(1.02,1),
                ncol=2, fontsize=8, frameon=False
            )
            plt.tight_layout(rect=[0,0,0.85,1])

            out = os.path.join(self.input_dir, f'{level}_barplot.pdf')
            fig.savefig(out, dpi=450, format='pdf')
            plt.close(fig)
            self.log(f"Saved {out}")

    def run_all(self):
        self.log("=== run_all start ===")
        for fp in self.file_paths:
            df = pd.read_csv(fp)
            self.filter_and_truncate(df, fp)
        self.compute_unclassified_stats()
        self.plot_well_classified()
        self.plot_taxa_retained()
        self.plot_cumulative_barplots()
        self.log("=== run_all complete ===")

if __name__ == '__main__':
    mt = MetabarcodingTaxonomy(
        input_dir='.',
        level_pattern='level-*.csv',
        level_labels=[
            'Kingdom', 'Phylum', 'Class',
            'Order', 'Family', 'Genus', 'Species'
        ]
    )
    mt.run_all()