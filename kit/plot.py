from numpy import arange
import matplotlib.pyplot as plt

class Plot:
    def __init__(self):
        self._data = []
        self._label = None
        self._fig = plt.figure()
        self._host = self._fig.add_subplot(111)
        self._labelFlag, self._xlabelFlag, self._ylabelFlag = False, False, False
        self._xticksFlag, self._yticksFlag = False, False
        plt.subplots_adjust(right=0.95, top=0.95)
    @property
    def data(self):
        return self._data
    @property
    def line_label(self):
        return self._labelFlag
    @line_label.setter
    def line_label(self, label):
        assert len(label) == len(self._data), "Label length does not match data length."
        self._label = label
        self._labelFlag = True
    @property
    def xlabel(self):
        return self._xlabelFlag
    @xlabel.setter
    def xlabel(self, xlabel):
        self._host.set_xlabel(xlabel, {"fontsize": 14})
        self._xlabelFlag = True
    @property
    def ylabel(self):
        return self._ylabelFlag
    @ylabel.setter
    def ylabel(self, ylabel):
        self._host.set_ylabel(ylabel, {"fontsize": 14})
        self._ylabelFlag = True
    @property
    def xticks(self):
        return self._xticksFlag
    @xticks.setter
    def xticks(self, xticks):
        data = []
        for iterator in range(len(self._data)):
            data += self._data[iterator][0]
        rng = max(data) - min(data)
        plt.xticks(arange(min(data), max(data)+rng/100, rng/(len(xticks)-1)), xticks, fontsize=12)
        self._xticksFlag = True
    @property
    def yticks(self):
        return self._yticksFlag
    @yticks.setter
    def yticks(self, yticks):
        data = []
        for iterator in range(len(self._data)):
            data += self._data[iterator][1]
        rng = max(data) - min(data)
        self._yticksFlag = True
        self._host.yticks(arange(min(data), max(data)+rng/100, rng/len(yticks)), yticks, fontsize=12)
    def auto_xticks(self):
        from numpy import arange, ceil
        x_len = self._data[0][0]
        if len(x_len) < 1001:
            xticks = None
        elif len(x_len) > 1000 and len(x_len) < 5001:
            xticks = list(arange(0, ceil(len(x_len/1000))))
        elif len(x_len) > 5000 and len(x_len) < 10001:
            xticks = list(arange(0, ceil(len(x_len)/1000), 2))
        elif len(x_len) > 10000 and len(x_len) < 25001:
            xticks = list(arange(0, ceil(len(x_len)/1000), 5))
        elif len(x_len) > 25000 and len(x_len) < 80001:
            xticks = list(arange(0, ceil(len(x_len)/1000), 10))
        elif len(x_len) > 80000 and len(x_len) < 200001:
            xticks = list(arange(0, ceil(len(x_len)/1000), 25))
        elif len(x_len) > 200000 and len(x_len) < 1000001:
            xticks = list(arange(0, ceil(len(x_len)/1000), 100))
        elif len(x_len) > 1000000:
            xticks = list(arange(0, ceil(len(x_len)/1000000)))
        if xticks is not None:
            self.xticks(xticks)
    def append(self, x, y):
        self._data.append((x, y))
    def plot(self, legend=True, palette=None):
        from seaborn import color_palette
        if palette is None:
            palette = color_palette("Set2", 8)
            if len(self.data) > 8:
                palette += color_palette("husl", len(self.data)-8)
        if len(self._data) < 6:
            fontsize = 14
        elif len(self._data) < 9:
            fontsize = 12
        elif len(self._data) < 11:
            fontsize = 10
        else:
            fontsize = 8
        if legend:
            for idx, data in enumerate(self._data):
                if self._label is None:
                    self._host.plot(data[0], data[1], color=palette[idx], label=str(idx), linewidth=2.2)
                else:
                    self._host.plot(data[0], data[1], color=palette[idx], label=self._label[idx], linewidth=2.2)
            self._host.legend(loc="best", fontsize=fontsize)
        if not self._xticksFlag:
            self._host.tick_params(axis='x', labelsize=12)
        if not self._yticksFlag:
            self._host.tick_params(axis='y', labelsize=12)
        for idx, data in enumerate(self._data):
            self._host.plot(data[0], data[1], color=palette[idx], linewidth=2.2)
        tmp_x, tmp_y = [], []
        for i in range(len(self._data)):
            tmp_x += list(self._data[i][0])
            tmp_y += list(self._data[i][1])
        buffer = (max(tmp_y) - min(tmp_y)) * 0.05
        self._host.set_xlim(min(tmp_x), max(tmp_x))
        self._host.set_ylim(min(tmp_y)-buffer, max(tmp_y)+buffer)
    def save(self, name):
        from kit.fundamental import Args
        from os import getcwd
        if ".png" not in name:
            name += ".png"
        Args.same_name(getcwd(), name)
        self._fig.savefig(name, dpi=600, format="png")

class TwinPlot(Plot):
    def __init__(self):
        Plot.__init__(self)
        self._left_data, self._right_data = [], []
        self._fig = plt.figure()
        self._host = self._fig.add_subplot(111)
        self._twin = self._host.twinx()
        plt.subplots_adjust(top=0.95)
    @property
    def data(self):
        return self._data
    @property
    def xlabel(self):
        return self._xlabelFlag
    @xlabel.setter
    def xlabel(self, xlabel):
        self._xlabelFlag = True
        self._host.set_xlabel(xlabel, {"fontsize": 14})
    @property
    def left_ylabel(self):
        return self._ylabelFlag % 2
    @left_ylabel.setter
    def left_ylabel(self, ylabel):
        self._ylabelFlag += 2
        self._host.set_ylabel(ylabel, {"fontsize": 14})
    @property
    def right_ylabel(self):
        return self._ylabelFlag // 2
    @right_ylabel.setter
    def right_ylabel(self, ylabel):
        self._ylabelFlag += 1
        self._twin.set_ylabel(ylabel, {"fontsize": 14})
    @property
    def left_yticks(self):
        return self._yticksFlag % 2
    @left_yticks.setter
    def left_yticks(self, yticks):
        data = []
        for iterator in range(len(self._left_data)):
            data += self._data[iterator][1]
        rng = max(data) - min(data)
        self._yticksFlag += 2
        self._host.yticks(arange(min(data), max(data)+rng/100, rng/len(yticks)), yticks, fontsize=14)
    @property
    def right_yticks(self):
        return self._yticksFlag // 2
    @right_yticks.setter
    def right_yticks(self, yticks):
        data = []
        for iterator in range(len(self._right_data)):
            data += self._data[iterator][1]
        rng = max(data) - min(data)
        self._yticksFlag += 1
        self._twin.yticks(arange(min(data), max(data)+rng/100, rng/len(yticks)), yticks, fontsize=14)
    def append(self, x, y, axis="left"):
        if axis == "left":
            self._left_data.append((x, y))
            self._data += self._left_data
        elif axis == "right":
            self._right_data.append((x, y))
            self._data += self._right_data
    def plot(self, palette=None, legend=True):
        from seaborn import color_palette
        if palette is None:
            palette = color_palette("Set2", 8)
            if len(self.data) > 8:
                palette += color_palette("husl", len(self.data)-8)
        if len(self._left_data) < 6:
            left_fontsize = 14
        elif len(self._left_data) < 9:
            left_fontsize = 12
        elif len(self._left_data) < 11:
            left_fontsize = 10
        else:
            left_fontsize = 8
        if len(self._right_data) < 6:
            right_fontsize = 14
        elif len(self._right_data) < 9:
            right_fontsize = 12
        elif len(self._right_data) < 11:
            right_fontsize = 10
        else:
            right_fontsize = 8
        if not self._xticksFlag:
            self._host.tick_params(axis='x', labelsize=12)
            self._twin.tick_params(axis='x', labelsize=12)
        if not self._yticksFlag:
            self._host.tick_params(axis='y', labelsize=12)
            self._twin.tick_params(axis='y', labelsize=12)
        if legend:
            for idx, data in enumerate(self._left_data):
                if len(self._label) == 0:
                    self._host.plot(data[0], data[1], color=palette[idx], label=str(idx), linewidth=2.2)
                else:
                    self._host.plot(data[0], data[1], color=palette[idx], label=self._label[idx], linewidth=2.2)
            for idx, data in enumerate(self._right_data, start=len(self._left_data)):
                if len(self._label) == 0:
                    self._twin.plot(data[0], data[1], color=palette[idx], label=str(idx), linewidth=2.2)
                else:
                    self._twin.plot(data[0], data[1], color=palette[idx], label=self._label[idx-len(self._left_data)+1], linewidth=2.2)
            self._host.legend(loc="upper left", fontsize=left_fontsize)
            self._twin.legend(loc="upper right", fontsize=right_fontsize)
        else:
            for idx, data in enumerate(self._left_data):
                self._host.plot(data[0], data[1], color=palette[idx], linewidth=2.2)
            for idx, data in enumerate(self._right_data, start=len(self._left_data)):
                self._twin.plot(data[0], data[1], color=palette[idx], linewidth=2.2)
        tmp_x, tmp_y_left, tmp_y_right = [], [], []
        for i in range(len(self._left_data)):
            tmp_x += list(self._data[i][0])
            tmp_y_left += list(self._left_data[i][1])
            tmp_y_right += list(self._right_data[i][1])
        left_buffer, right_buffer = (max(tmp_y_left) - min(tmp_y_left)) * 0.05, (max(tmp_y_right) - min(tmp_y_right)) * 0.05
        self._host.set_xlim(min(tmp_x), max(tmp_x))
        self._twin.set_xlim(min(tmp_x), max(tmp_x))
        self._host.set_ylim(min(tmp_y_left)-left_buffer, max(tmp_y_left)+left_buffer)
        self._twin.set_ylim(min(tmp_y_right)-right_buffer, max(tmp_y_right)+right_buffer)

def plot(fig, file_name, ticks=False, line_label=True, xlabel=True, ylabel=True, palette=None):
    if ticks:
        from numpy import arange, ceil
        while(1):
            tmp = input("Adjust X axis ticks, or auto mode (y/n/a) [a]: ")
            if tmp == '' or tmp.lower() == 'a':
                x_len = fig.data[0][0]
                if len(x_len) > 1000 and len(x_len) < 5001:
                    fig.xticks = list(arange(0, ceil(len(x_len/1000))))
                elif len(x_len) > 5000 and len(x_len) < 10001:
                    fig.xticks = list(arange(0, ceil(len(x_len)/1000), 2))
                elif len(x_len) > 10000 and len(x_len) < 25001:
                    fig.xticks = list(arange(0, ceil(len(x_len)/1000), 5))
                elif len(x_len) > 25000 and len(x_len) < 80001:
                    fig.xticks = list(arange(0, ceil(len(x_len)/1000), 10))
                elif len(x_len) > 80000 and len(x_len) < 200001:
                    fig.xticks = list(arange(0, ceil(len(x_len)/1000), 25))
                elif len(x_len) > 200000 and len(x_len) < 1000001:
                    fig.xticks = list(arange(0, ceil(len(x_len)/1000), 100))
                elif len(x_len) > 1000000:
                    fig.xticks = list(arange(0, ceil(len(x_len)/1000000)))
                break
            elif tmp.lower() == 'n':
                break
            else:
                print("Warning: Input error.")
    if not line_label:
        nameFlag = False
    if line_label and not fig.line_label:
        while(1):
            tmp = input("Name the lines (y/n) [n]: ")
            if tmp == '':
                nameFlag = False
                break
            elif tmp.lower() == 'y':
                nameFlag = True
                break
            elif tmp.lower() == 'n':
                nameFlag = False
                break
            else:
                print("Warning: Input error.")
        if nameFlag:
            label = []
            if isinstance(fig, Plot):
                for i in range(len(fig.data)):
                    line = input(f"Line {i+1}: ")
                    label.append(line)
            elif isinstance(fig, TwinPlot):
                for i in range(len(fig.data[0])):
                    line = input(f"Left line {i+1}: ")
                    label.append(line)
                    line = input(f"Right line {i+1}: ")
                    label.append(line)
            fig.line_label = label
    if xlabel and not fig.xlabel:
        while(1):
            tmp = input("Name X axis (y/n) [n]: ")
            if tmp == '':
                break
            elif tmp.lower() == 'y':
                fig.xlabel = input("X axis: ")
                break
            elif tmp.lower() == 'n':
                break
            else:
                print("Warning: Input error.")
    if isinstance(fig, Plot):
        if not fig.ylabel:
            while(1):
                tmp = input("Name Y axis (y/n) [n]: ")
                if tmp == '':
                    break
                elif tmp.lower() == 'y':
                    fig.ylabel = input("Y axis: ")
                    break
                elif tmp.lower() == 'n':
                    break
                else:
                    print("Warning: Input error.")
    elif ylabel and isinstance(fig, TwinPlot):
        if not fig.left_ylabel:
            while(1):
                tmp = input("Name left Y axis (y/n) [n]: ")
                if tmp == '':
                    break
                elif tmp.lower() == 'y':
                    fig.left_ylabel = input("Left Y axis: ")
                    break
                elif tmp.lower() == 'n':
                    break
                else:
                    print("Warning: Input error.")
        elif ylabel and not fig.right_ylabel:
            while(1):
                tmp = input("Name right Y axis (y/n) [n]: ")
                if tmp == '':
                    break
                elif tmp.lower() == 'y':
                    fig.right_ylabel = input("Right Y axis: ")
                    break
                elif tmp.lower() == 'n':
                    break
                else:
                    print("Warning: Input error.")
    if palette is None:
        while(1):
            tmp = input("Adjust color (y/n) [n]: ")
            if tmp == '':
                palette = None
                break
            elif tmp.lower() == 'y':
                palette = []
                for i in range(len(fig.data)):
                    tmp = input(f"Color {i+1}: ")
                    palette.append(tmp)
                break
            elif tmp.lower() == 'n':
                palette = None
                break
            else:
                print("Warning: Input error.")
    fig.plot(palette=palette, legend=nameFlag)
    fig.save(file_name)
