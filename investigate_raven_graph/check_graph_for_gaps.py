import argparse
import dgl

def interval_union(path):
    graph = dgl.load_graphs(path)[0][0]
    intervals = []
    for strand, start, end in zip(graph.ndata['read_strand'], graph.ndata['read_start'], graph.ndata['read_end']):
        if strand.item() == 1:
            intervals.append([start.item(), end.item()])
    intervals.sort(key=lambda x: x[0])
    result = [intervals[0]]

    for interval in intervals[1:]:
        if interval[0] <= result[-1][1]:
            result[-1][1] = max(result[-1][1], interval[1])
        else:
            result.append(interval)

    return result

def run(args):
    # path = f'{root}/processed/{name}.dgl'
    intervals = interval_union(args.path)
    print(intervals)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--path', type=str, default='data/chr22_raven.dgl', help='Path to dgl file')


    args = parser.parse_args()
    run(args)