Installation
============

Pip-Install
-------
It is recommended to use a python environment.

.. code-block:: console

    # Create a python environment (optional) or use an existing one
    python -m venv ~/venv/sasa
    # Activate your Python environment
    source ~/venv/sasa/bin/activate

From this projects root directory install the package with all dependencies

.. code-block:: console

    pip install .

Verification
------------

After installation, you can either verify that the SASA extension is loaded correctly (tp not confuse python with the local `sasa_ext` directory, execute this from any other than the source-file directory):

.. code-block:: console

    cd && python -c "import sasa_ext; print('✓ SASA extension loaded successfully')"

or/and execute the tests in the `tests` directory. This requires the installation of the test-requirements and the install in editable mode:

.. code-block:: console
    
    pip install -e ".[test]"
    pytest ./tests
