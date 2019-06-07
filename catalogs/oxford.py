import zarr


def ag1000g_phase1_ar3_variation_main():
    """Open the Ag1000G phase 1 AR3 main variation callset.

    Returns
    -------
    zarr.Group

    """
    return zarr.open_consolidated(
        "/kwiat/vector/ag1000g/release/phase1.AR3/variation/main/zarr/ag1000g.phase1.ar3"
    )


def ag1000g_phase1_ar3_variation_main_pass():
    """Open the Ag1000G phase 1 AR3 main variation callset, pass variants subset.

    Returns
    -------
    zarr.Group

    """
    return zarr.open_consolidated(
        "/kwiat/vector/ag1000g/release/phase1.AR3/variation/main/zarr/ag1000g.phase1.ar3.pass"
    )


def ag1000g_phase1_ar3_variation_main_pass_biallelic():
    """Open the Ag1000G phase 1 AR3 main variation callset, pass biallelic
    variants subset.

    Returns
    -------
    zarr.Group

    """
    return zarr.open_consolidated(
        "/kwiat/vector/ag1000g/release/phase1.AR3/variation/main/zarr/ag1000g.phase1.ar3.pass.biallelic"
    )


def ag1000g_phase2_ar1_variation_main():
    """Open the Ag1000G phase 2 AR1 main variation callset.

    Returns
    -------
    zarr.Group

    """
    return zarr.open_consolidated(
        "/kwiat/vector/ag1000g/release/phase2.AR1/variation/main/zarr/all/ag1000g.phase2.ar1"
    )


def ag1000g_phase2_ar1_variation_main_pass():
    """Open the Ag1000G phase 2 AR1 main variation callset, pass variants subset.

    Returns
    -------
    zarr.Group

    """
    return zarr.open_consolidated(
        "/kwiat/vector/ag1000g/release/phase2.AR1/variation/main/zarr/pass/ag1000g.phase2.ar1.pass"
    )


def ag1000g_phase2_ar1_variation_main_pass_biallelic():
    """Open the Ag1000G phase 2 AR1 main variation callset, pass biallelic
    variants subset.

    Returns
    -------
    zarr.Group

    """
    return zarr.open_consolidated(
        "/kwiat/vector/ag1000g/release/phase2.AR1/variation/main/zarr/biallelic/ag1000g.phase2.ar1.pass.biallelic"
    )
