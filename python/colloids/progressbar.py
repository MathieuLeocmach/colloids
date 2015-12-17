#Code from https://github.com/ipython/ipython/issues/1527/

import sys, time
try:
    from IPython.core.display import clear_output
    have_ipython = True
    try:
        from ipykernel.zmqshell import ZMQInteractiveShell
        from ipywidgets import FloatProgress
        from IPython.display import display
        in_notebook = isinstance(get_ipython(), ZMQInteractiveShell)
    except ImportError:
        in_notebook = False
except ImportError:
    have_ipython = False

class ProgressBar:
    """Displays the progression of N steps known in advance (i.e. a loop of N iterations).   Example of use:
    pro = ProgressBar(5000)
    for i in range(5000):
        #do something
        pro.animate(i)
    """
    def __init__(self, iterations):
        self.iterations = iterations
        self.prog_bar = '[]'
        self.fill_char = '*'
        self.width = 40
        self.__update_amount(0)
        if have_ipython:
            if(in_notebook):
                self.animate = self.animate_notebook
                self.fp = FloatProgress(min=0, max=self.iterations)
                display(self.fp)
            else:
                self.animate = self.animate_ipython
        else:
            self.animate = self.animate_noipython

    def animate_ipython(self, iter):
        try:
            clear_output()
        except Exception:
            # terminal IPython has no clear_output
            pass
        print(self, end='\r')
        sys.stdout.flush()
        self.update_iteration(iter + 1)
        
    def animate_notebook(self, i):
        self.fp.value = i

    def update_iteration(self, elapsed_iter):
        self.__update_amount((elapsed_iter / float(self.iterations)) * 100.0)
        self.prog_bar += '  %d of %s complete' % (elapsed_iter, self.iterations)

    def __update_amount(self, new_amount):
        percent_done = int(round((new_amount / 100.0) * 100.0))
        all_full = self.width - 2
        num_hashes = int(round((percent_done / 100.0) * all_full))
        self.prog_bar = '[' + self.fill_char * num_hashes + ' ' * (all_full - num_hashes) + ']'
        pct_place = (len(self.prog_bar) // 2) - len(str(percent_done))
        pct_string = '%d%%' % percent_done
        self.prog_bar = self.prog_bar[0:pct_place] + \
            (pct_string + self.prog_bar[pct_place + len(pct_string):])

    def __str__(self):
        return str(self.prog_bar)
