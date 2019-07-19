import concurrent.futures
import os

MAX_PROCESS_NUMBER = os.cpu_count()
MAX_JOB_NUMBER = MAX_PROCESS_NUMBER + 1


def launch_multiprocess(runs, fn):
    """
    Spawn multiprocesses avoiding memory problems.
    The 'runs' parameter must be a list of 'cmd' s.t., subprocess.Popen(cmd).wait() works.
    :param runs: a list of cmd (parameters, executable, etc).
    :param fn: the function to be parallelized
    :return:
    """
    done = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=MAX_PROCESS_NUMBER) as executor:
        jobs = {}
        runs_left = len(runs)
        runs_iter = iter(runs)

        while runs_left:
            for run in runs_iter:
                job = executor.submit(fn, run)
                jobs[job] = run
                if len(jobs) > MAX_JOB_NUMBER:
                   break

            for job in concurrent.futures.as_completed(jobs):
                runs_left -= 1
                result = job.result()
                run = jobs[job]
                del jobs[job]
                done.append(result)
                break
    return done
