from common.common import *

# -----------------------------------------------------------------------------
def test_safe_figure_save():
# -----------------------------------------------------------------------------
    import os
    import matplotlib.pylab as pyl
    import numpy as np

    path = "/tmp/"
    file_name = "file_name.png"

    pyl.plot(np.linspace(0,10,6))
    mkdir_p(path)
    
    open(path + file_name, "a").close() 
    
    if os.path.exists("/tmp/file_name_0.png"):
        os.remove("/tmp/file_name_0.png")

    save_figure_safe(path+file_name)

    assert os.path.exists("/tmp/file_name_0.png")
   
    # Clean-up after the test
    os.remove("/tmp/file_name_0.png")
    os.remove("/tmp/file_name.png")
    
# -----------------------------------------------------------------------------
def test_OU_process():
# -----------------------------------------------------------------------------
    """
    Trying out OU_process method by comparing it to a method I found online.
    Since the one online uses different random generator, it will be impossible
    to make a direct comparison ... Not sure how to handle this.
    """
    pass
    