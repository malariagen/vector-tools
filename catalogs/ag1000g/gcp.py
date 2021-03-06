import zarr
from pathlib import Path
from cached_property import cached_property
import gcsfs


GCP_PROJECT = 'malariagen-jupyterhub'
AG1000G_RELEASE_DIR = Path("ag1000g-release")
PHASE1_AR3_DIR = AG1000G_RELEASE_DIR / 'phase1.AR3'
PHASE1_AR31_DIR = AG1000G_RELEASE_DIR / 'phase1.AR3.1'
PHASE2_AR1_DIR = AG1000G_RELEASE_DIR / 'phase2.AR1'


class Phase1AR3(object):

    def __init__(self):
        self.fs = gcsfs.GCSFileSystem(project=GCP_PROJECT, token='anon',
                                      access='read_only')

    @cached_property
    def variation_main(self):
        path = PHASE1_AR3_DIR / 'variation/main/zarr/ag1000g.phase1.ar3'
        store = gcsfs.GCSMap(str(path), gcs=self.fs, check=False, create=False)
        return zarr.open_consolidated(store)

    @cached_property
    def variation_main_pass(self):
        path = PHASE1_AR3_DIR / 'variation/main/zarr/ag1000g.phase1.ar3.pass'
        store = gcsfs.GCSMap(str(path), gcs=self.fs, check=False, create=False)
        return zarr.open_consolidated(store)

    @cached_property
    def variation_main_pass_biallelic(self):
        path = PHASE1_AR3_DIR / 'variation/main/zarr/ag1000g.phase1.ar3.pass.biallelic'
        store = gcsfs.GCSMap(str(path), gcs=self.fs, check=False, create=False)
        return zarr.open_consolidated(store)


class Phase1AR31(object):

    def __init__(self):
        self.fs = gcsfs.GCSFileSystem(project=GCP_PROJECT, token='anon',
                                      access='read_only')

    @cached_property
    def haplotypes_main(self):
        path = PHASE1_AR31_DIR / 'haplotypes/main/zarr/ag1000g.phase1.ar3.1.haplotypes'
        store = gcsfs.GCSMap(str(path), gcs=self.fs, check=False, create=False)
        return zarr.open_consolidated(store)


class Phase2AR1(object):

    def __init__(self):
        self.fs = gcsfs.GCSFileSystem(project=GCP_PROJECT, token='anon',
                                      access='read_only')

    @cached_property
    def variation_main(self):
        path = PHASE2_AR1_DIR / 'variation/main/zarr/all/ag1000g.phase2.ar1'
        store = gcsfs.GCSMap(str(path), gcs=self.fs, check=False, create=False)
        return zarr.open_consolidated(store)

    @cached_property
    def variation_main_pass(self):
        path = PHASE2_AR1_DIR / 'variation/main/zarr/pass/ag1000g.phase2.ar1.pass'
        store = gcsfs.GCSMap(str(path), gcs=self.fs, check=False, create=False)
        return zarr.open_consolidated(store)

    @cached_property
    def variation_main_pass_biallelic(self):
        path = PHASE2_AR1_DIR / 'variation/main/zarr/biallelic/ag1000g.phase2.ar1.pass.biallelic'
        store = gcsfs.GCSMap(str(path), gcs=self.fs, check=False, create=False)
        return zarr.open_consolidated(store)

    @cached_property
    def haplotypes_main(self):
        path = PHASE2_AR1_DIR / 'haplotypes/main/zarr/ag1000g.phase2.ar1.haplotypes'
        store = gcsfs.GCSMap(str(path), gcs=self.fs, check=False, create=False)
        return zarr.open_consolidated(store)
