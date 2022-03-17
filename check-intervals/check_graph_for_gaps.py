import argparse
import dgl
import interval

def create_inteval_list(path):
    graph = dgl.load_graphs(path)[0][0]
    intervals = []
    for strand, start, end in zip(graph.ndata['read_strand'], graph.ndata['read_start'], graph.ndata['read_end']):
        if strand.item() == 1:
            intervals.append([start.item(), end.item()])
    return intervals

def run(args):
    # path = f'{root}/processed/{name}.dgl'
    interval_list = create_inteval_list(args.path)
    intervals = interval.interval_union(interval_list)
    print(intervals)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--path', type=str, default='data/chr22_raven.dgl', help='Path to dgl file')

    args = parser.parse_args()
    run(args)