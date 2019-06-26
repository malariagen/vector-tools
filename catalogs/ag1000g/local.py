import zarr
from pathlib import Path
from cached_property import cached_property


# default base directories (can be overridden)
AG1000G_RELEASE_DIR = Path("/kwiat/vector/ag1000g/release")
PHASE1_AR3_DIR = AG1000G_RELEASE_DIR / 'phase1.AR3'
PHASE1_AR31_DIR = AG1000G_RELEASE_DIR / 'phase1.AR3.1'
PHASE2_AR1_DIR = AG1000G_RELEASE_DIR / 'phase2.AR1'


class Phase1AR3(object):

    def __init__(self, base_dir=PHASE1_AR3_DIR):
        self.base_dir = Path(base_dir)

    @cached_property
    def variation_main(self):
        path = self.base_dir / 'variation/main/zarr/ag1000g.phase1.ar3'
        return zarr.open_consolidated(str(path))

    @cached_property
    def variation_main_pass(self):
        path = self.base_dir / 'variation/main/zarr/ag1000g.phase1.ar3.pass'
        return zarr.open_consolidated(str(path))

    @cached_property
    def variation_main_pass_biallelic(self):
        path = self.base_dir / 'variation/main/zarr/ag1000g.phase1.ar3.pass.biallelic'
        return zarr.open_consolidated(str(path))


class Phase1AR31(object):

    def __init__(self, base_dir=PHASE1_AR31_DIR):
        self.base_dir = Path(base_dir)

    @cached_property
    def haplotypes_main(self):
        path = self.base_dir / 'haplotypes/main/zarr/ag1000g.phase1.ar3.1.haplotypes'
        return zarr.open_consolidated(str(path))

    
class Phase2AR1(object):

    def __init__(self, base_dir=PHASE2_AR1_DIR):
        self.base_dir = Path(base_dir)

    @cached_property
    def variation_main(self):
        path = self.base_dir / 'variation/main/zarr/all/ag1000g.phase2.ar1'
        return zarr.open_consolidated(str(path))

    @cached_property
    def variation_main_pass(self):
        path = self.base_dir / 'variation/main/zarr/pass/ag1000g.phase2.ar1.pass'
        return zarr.open_consolidated(str(path))

    @cached_property
    def variation_main_pass_biallelic(self):
        path = self.base_dir / 'variation/main/zarr/biallelic/ag1000g.phase2.ar1.pass.biallelic'
        return zarr.open_consolidated(str(path))

    @cached_property
    def haplotypes_main(self):
        path = self.base_dir / 'haplotypes/main/zarr/ag1000g.phase2.ar1.haplotypes'
        return zarr.open_consolidated(str(path))
