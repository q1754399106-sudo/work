"""Read-only snapshot of the running one-shot research pipeline."""
from pathlib import Path
import csv
import datetime as dt
import json

ROOT=Path(__file__).resolve().parent


def read(path):
    try:
        return json.loads(path.read_text(encoding='utf-8-sig'))
    except (OSError,ValueError):
        return {}


def snapshot():
    pipeline=read(ROOT/'experiments/test_v2_PIPELINE.json')
    state=read(ROOT/'experiments/test_v2/matlab/run_status.json')
    last={}
    path=ROOT/'experiments/test_v2/matlab/metrics_partial.csv'
    try:
        with path.open(encoding='utf-8-sig',newline='') as stream:
            for row in csv.DictReader(stream):
                last=row
    except (OSError,csv.Error):
        pass
    return dict(utc=dt.datetime.now(dt.timezone.utc).isoformat(),pipeline=pipeline.get('status'),
                stage=pipeline.get('active_stage'),completed_observations=state.get('observations'),
                completed_rows=state.get('rows'),failed_rows=state.get('failed_rows'),
                elapsed_seconds=state.get('elapsed_seconds'),
                latest_case={k:last.get(k) for k in ['Model','Suite','Scenario','CalibrationSeed','Seed']})


if __name__=='__main__':
    print(json.dumps(snapshot(),ensure_ascii=True))
