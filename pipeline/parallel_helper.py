from functools import partial
from multiprocessing import Pool
from tqdm import tqdm
import sys

def run_parallel(func, args, n_processes = 1, chunk_callback=None, return_results = False, **kwargs):
    
    p = Pool(n_processes)
    results = []

    if len(kwargs) > 0:
        func = partial(func, **kwargs)

    with tqdm(total = len(args)) as pbar:
        for res in tqdm(p.imap_unordered(func, args, chunksize=n_processes)):
            pbar.update()
            if type(res) == bool or res == None:
                continue

            if chunk_callback is not None:
                try:
                    chunk_callback(res, **kwargs)
                except Exception as e:
                    print("Error during insertion to db", sys.exc_info()[0], e)

            if return_results:
                results.append(res)

    pbar.close()
    p.close()
    p.join()
    
    return results if return_results else True