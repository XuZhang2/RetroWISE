def run_imap_mp(func, argument_list, num_processes='', is_tqdm=True):
    '''
    多进程与进度条结合

    param:
    ------
    func:function
        函数
    argument_list:list
        参数列表
    num_processes:int
        进程数，不填默认为总核心-3
    is_tqdm:bool
        是否展示进度条，默认展示
    '''
    result_list_tqdm = []
    try:
        import multiprocessing
        if num_processes == '':
            num_processes = multiprocessing.cpu_count()-3
        pool = multiprocessing.Pool(processes=num_processes)
        if is_tqdm:
            from tqdm import tqdm
            for result in tqdm(pool.imap(func=func, iterable=argument_list), total=len(argument_list)):
                result_list_tqdm.append(result)
        else:
            for result in pool.imap(func=func, iterable=argument_list):
                result_list_tqdm.append(result)
        pool.close()
        pool.join()
    except:
        result_list_tqdm = list(map(func,argument_list))
    return result_list_tqdm
