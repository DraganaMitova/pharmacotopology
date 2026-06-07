from __future__ import annotations

import urllib.error

from pharmacotopology.folding_msa_free_internet_predictors import post_esm_atlas_fold_sequence


def test_esm_atlas_adapter_fails_closed_without_persisting_sequence(monkeypatch, tmp_path):
    def _fake_urlopen(*args, **kwargs):  # noqa: ANN001
        raise urllib.error.URLError("offline test")

    monkeypatch.setattr("urllib.request.urlopen", _fake_urlopen)
    run = post_esm_atlas_fold_sequence(
        sequence="ACDEFGHIKLMNPQRSTVWY",
        output_pdb_path=tmp_path / "prediction.pdb",
        timeout_seconds=1,
        minimum_ca_fraction=0.8,
    )

    assert run.attempted is True
    assert run.success is False
    assert run.status == "request_failed"
    assert run.raw_sequence_persisted is False
    assert run.alphafold_endpoint_used is False
    assert not (tmp_path / "prediction.pdb").exists()
