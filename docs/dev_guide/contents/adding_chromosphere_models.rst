.. _adding_chromosphere_models:

===========================================
Adding a New Chromosphere Model to SunPy
===========================================

This guide provides instructions on how to add a new 1D chromosphere model to `sunpy.sun.models`.
SunPy provides pre-defined atmospheric models, such as `Avrett & Loeser (2008)`. These models are stored as structured data files and dynamically loaded into the `sunpy.sun.models` module.

----------------------------
Preparing the Model Data
----------------------------

Ensure your model data is in a structured format, preferably a QTable with the relevant parameters for model.
Save the data as an ECSV  file and place it in the `data/` directory within `sunpy/sun/`.

----------------------------
Implementing the Model
----------------------------

Modify `sunpy/sun/models.py` to include the new model:

.. code-block:: python

    import pathlib

    _MODEL_DATA_DIR = pathlib.Path(__file__).parent.absolute() / "data"
    _MODELS = {
        "chromosphere_avrett_loeser_2008": _MODEL_DATA_DIR / "chromosphere_avrett_loeser_2008_model.ecsv",
        "chromosphere_new_model": _MODEL_DATA_DIR / "chromosphere_new_model.ecsv",  # Add your model here
    }

This registers the model inside `sunpy.sun.models`.

----------------------------
Testing
----------------------------

- Ensure all tests pass by running `pytest`.
