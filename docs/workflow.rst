.. include:: links.rst

.. _Workflow:

HyDe Workflow
=============

.. code:: py

  import phyde as hd

  # Analyze the data in the `test/` folder within the main HyDe directory
  res, boot = hd.run_hyde("data.txt", "map.txt", "out", 16, 4, 50000, bootReps=100)

  # Filter the main results in the `res` variable to only include significant
  # results and gamma values between 0 and 1.
  # We can do this using a dicionary comprehension.
  filt_res = {k: v for k,v in res.items if v['Pvalue'] < 0.05 and v['Gamma'] > 0 and v['Gamma'] < 1}
