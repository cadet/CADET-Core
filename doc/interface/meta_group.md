(ffmeta)=

# Meta Group

`FILE_FORMAT`

> Version of the file format (defaults to 040000 = 4.0.0 if omitted) with two digits per part (Major.Minor.Patch)
>
> |   |   |
> | --- | --- |
> | **In/out:** In |  |

`CADET_VERSION`

> Version of the executed `CADET` simulator
>
> |   |   |
> | --- | --- |
> | **In/out:** Out |  |

`CADET_COMMIT`

> Git commit SHA1 from which the `CADET` simulator was built
>
> |   |   |
> | --- | --- |
> | **In/out:** Out |  |

`CADET_BRANCH`

> Git branch from which the `CADET` simulator was built
>
> |   |   |
> | --- | --- |
> | **In/out:** Out |  |

`TIME_SIM`

> Time that the time integration took (excluding any preparations and postprocessing)
>
> **Unit:** $\mathrm{s}$
>
> |   |   |
> | --- | --- |
> | **In/out:** Out |  |
