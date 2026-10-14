.. _FFMeta:

Meta Group
==========

``FILE_FORMAT``

   Version of the file format (defaults to 040000 = 4.0.0 if omitted) with two digits per part (Major.Minor.Patch)
   
   ================  =========================
   **In/out:** In    **Type:** int
   ================  =========================
   
``CADET_VERSION``

   Version of the executed :math:`\texttt{CADET}` simulator
   
   ================  =========================
   **In/out:** Out   **Type:** string
   ================  =========================
   
``CADET_COMMIT``

   Git commit SHA1 from which the :math:`\texttt{CADET}` simulator was built
   
   ================  =========================
   **In/out:** Out   **Type:** string
   ================  =========================
   
``CADET_BRANCH``

   Git branch from which the :math:`\texttt{CADET}` simulator was built
   
   ================  =========================
   **In/out:** Out   **Type:** string
   ================  =========================
   
``TIME_SIM``

   Time that the time integration took (excluding any preparations and postprocessing)

   **Unit:** :math:`\mathrm{s}`
   
   ================  =========================
   **In/out:** Out   **Type:** double
   ================  =========================

Group /meta/idas
----------------

Some of the optional outputs, especially the various counters, can be very useful in determining how successful the IDA solver is in doing its job.
For example, the counters nsteps and nrevals provide a rough measure of the overall cost of a given run, and can be compared among runs with differing input options to suggest which set of options is most efficient.
The ratio nniters/nsteps measures the performance of the nonlinear solver in solving the nonlinear systems at each time step; typical values for this range from 1.1 to 1.8.
The ratio njevals/nniters (in the case of a matrix-based linear solver), and the ratio npevals/nniters (in the case of an iterative linear solver) measure the overall degree of nonlinearity in these systems, and also the quality of the approximate Jacobian or preconditioner being used.
Thus, for example, njevals/nniters can indicate if a user-supplied Jacobian is inaccurate, if this ratio is larger than for the case of the corresponding internal Jacobian.
The ratio nliters/nniters measures the performance of the Krylov iterative linear solver, and thus (indirectly) the quality of the preconditioner.
In the following, `NSOLUTIONTIMES` denotes the number of time points at which the solution is reported, see `USER_SOLUTION_TIMES`.

``IDAS_NTIMESTEPS``

   Number of time steps taken by the IDAS time integrator

   ================  =========================
   **In/out:** Out   **Type:** int
   ================  =========================

``IDAS_NEXT_STEP_SIZE``

   Time integration step size to be attempted on the next internal step

   ================  =========================
   **In/out:** Out   **Type:** double
   ================  =========================

``IDAS_NEXT_TIME_POINT``

   Time integration next time point to be attempted on the next internal step

   ================  =========================
   **In/out:** Out   **Type:** double
   ================  =========================

``IDAS_BDF_ORDER``

    Integration method order at corresponding time point

   	================  ======================  ===============================
	**In/out:** Out   **Type:** int           **Length:** ``NSOLUTIONTIMES``
	================  ======================  ===============================

``IDAS_KUMULATIVE_NTIMESTEPS``

    Cumulative number of time steps taken by the IDAS at corresponding time point

   	================  ======================  ===============================
	**In/out:** Out   **Type:** int           **Length:** ``NSOLUTIONTIMES``
	================  ======================  ===============================

``IDAS_N_RESIDUAL_CALLS``

    Number of calls to the user's residual evaluation function at corresponding time point

   	================  ======================  ===============================
	**In/out:** Out   **Type:** int           **Length:** ``NSOLUTIONTIMES``
	================  ======================  ===============================

``IDAS_N_LINSOLVER_SETUP_CALLS``

    Cumulative number of calls made to the linear solver's setup function at corresponding time point

   	================  ======================  ===============================
	**In/out:** Out   **Type:** int           **Length:** ``NSOLUTIONTIMES``
	================  ======================  ===============================

``IDAS_N_LOCAL_ERROR_TEST_FAILURES``

    Cumulative number of local error test failures that have occurred at corresponding time point

   	================  ======================  ===============================
	**In/out:** Out   **Type:** int           **Length:** ``NSOLUTIONTIMES``
	================  ======================  ===============================

``IDAS_INTEGRATION_STEP_SIZE``

    Integration step size taken on the last internal step at corresponding time point

   	================  ======================  ===============================
	**In/out:** Out   **Type:** int           **Length:** ``NSOLUTIONTIMES``
	================  ======================  ===============================

``IDAS_FIRST_INTEGRATION_STEP_SIZE``

    Value of the integration step size used on the first step

   	================  ======================  ===============================
	**In/out:** Out   **Type:** int           **Length:** ``NSOLUTIONTIMES``
	================  ======================  ===============================
