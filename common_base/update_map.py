import sys

class UpdateMap():
    def __init__(self):
        self.map = {}

    def reply(self, worker_id, code):
        # Each worker must also reply
        # That it has updated successfully
        if code == 0:
            try:
                del self.map[worker_id]
            except KeyError:
                sys.stderr.write("Warning: not expecting worker %s to reply\n" % worker_id)
        else:
            sys.stderr.write("Warning: worker %s returned non-zero" % worker_id)
        
    def get(self, worker_id):
        # Each worker has to get mail
        if worker_id not in self.map:
            return None
        else:
            return self.map[worker_id].pop(0)
        print(self.map)

    def put(self, worker_id, update_str, *args):
        print("putting", worker_id, update_str)
        if worker_id in self.map:
            self.map[worker_id].append((update_str, *args))
        else:
            self.map[worker_id] = [(update_str, *args)]

    def is_empty(self):
        if len(self.map) == 0:
            return True
        return False

    def collect(self):
        while self.is_empty():
            pass
