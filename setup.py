from setuptools import setup, find_packages
  
setup(name='yiiplib',
      version="0.1.0",
      packages=find_packages(),
      package_data={"yiiplib": ["topology/5VRF_holo.pdb",
                                "topology/5VRF_of.pdb",
                                "topology/7KZX_apo.pdb",
                                "topology/7KZZ_holo.pdb"]},
      install_requires=["numpy",
                        "matplotlib",
                        "MDAnalysis",
                        "xarray",
                        "sklearn",
                        "seaborn",
                        "numkit",
                        "netcdf4"],
)