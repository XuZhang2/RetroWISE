def run_imap_mp(func, argument_list, num_processes='', is_tqdm=True):
    result_list_tqdm = []
    import multiprocessing
    if num_processes == '':
        num_processes = multiprocessing.cpu_count()-3
    with multiprocessing.Pool(processes=num_processes) as pool:  # Use context manager here
        if is_tqdm:
            from tqdm import tqdm
            for result in tqdm(pool.imap(func=func, iterable=argument_list), total=len(argument_list)):
                result_list_tqdm.append(result)
        else:
            for result in pool.imap(func=func, iterable=argument_list):
                result_list_tqdm.append(result)
            # pool.close()  # Not needed with the context manager
            # pool.join()  # Not needed with the context manager
    return result_list_tqdm

